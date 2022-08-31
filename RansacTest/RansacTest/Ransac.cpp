#include "pch.h"
#include "Ransac.h"

Ransac::Ransac()
{
	srand(getTickCount());
}

void Ransac::Clear()
{
//	m_score = 0;
	m_best.clear();
}

// 주성분 분석 사용모델
bool Ransac::getBestModel(const vector<Point>& total, const double percent, const double distance)
{
	const int max_pair = total.size() * percent;
	vector<int> sample_num;
	sample_num.reserve(max_pair * 2);

	// 이전과 중복 데이터 인지 확인한다. (떨림방지)
	if (total.size() == m_old.size())
	{
		if (total == m_old) return true;
	}

	// make random sample
	for (int i = 0; i < max_pair; i++)
	{
		int mod_num = total.size();

		int random1 = rand() % mod_num;
		int random2 = rand() % mod_num;

		if (random1 != random2)
		{
			sample_num.push_back(random1);
			sample_num.push_back(random2);
		}
		else i--;	// try again
	}

	// 
	int big_inliercount = 0;

	vector<Point> sample;
	sample.reserve(max_pair * 2);
	m_best.reserve(max_pair * 2);

	for (int i = 0; i < sample_num.size(); i += 2)
	{
		const int p1_idx = sample_num[i];
		const int p2_idx = sample_num[i + 1];

		Point p1 = total[p1_idx];
		Point p2 = total[p2_idx];

		// 두점의 직선을 구한다. ax+by+c=0
		double a = p1.y - p2.y;
		double b = p2.x - p1.x;
		double c = p1.x * p2.y - p2.x * p1.y;

		// Return distance between passed "point" and this line

		const double denominator = sqrt(a * a + b * b);	// 거리 계산에 쓰이는 분모

		sample.clear();
		int inlier_count = 0;

		// 모든 점을 대상으로 거리를 구한다.
		for (int k = 0; k < total.size(); k++)
		{
			Point p = total[k];
			double d = fabs(p.x * a + p.y * b + c) / denominator;

			if (d < distance)
			{
				inlier_count++;
				sample.push_back(p);
			}
		}

		// save best sample
		if (big_inliercount < inlier_count)
		{
			big_inliercount = inlier_count;
			m_best.swap(sample);	// fast m_best = sample;
		}

	}

	// PCA 모델 계산
	compute_model(m_best);

	// 라인 점수: 모든 점들 중에 inlier 비율
	m_score = (float)big_inliercount / total.size();
	m_old = total;
	return true;

}

// PCA opencv 계산
int Ransac::compute_model(const vector<Point>& pts)
{
	//Construct a buffer used by the pca analysis
	Mat data_pts = Mat(pts.size(), 2, CV_64F);
	for (int i = 0; i < data_pts.rows; i++)
	{
		data_pts.at<double>(i, 0) = pts[i].x;
		data_pts.at<double>(i, 1) = pts[i].y;
	}
	//Perform PCA analysis
	PCA pca_analysis(data_pts, Mat(), PCA::DATA_AS_ROW);
	//Store the center of the object
	Point cntr = Point(static_cast<int>(pca_analysis.mean.at<double>(0, 0)),
		static_cast<int>(pca_analysis.mean.at<double>(0, 1)));
	
	//Store the eigenvalues and eigenvectors
	vector<Point2d> eigen_vecs(2);
	vector<double> eigen_val(2);
	for (int i = 0; i < 2; i++)
	{
		eigen_vecs[i] = Point2d(pca_analysis.eigenvectors.at<double>(i, 0),
			pca_analysis.eigenvectors.at<double>(i, 1));
		eigen_val[i] = pca_analysis.eigenvalues.at<double>(i);
	}



	// 제1성분 벡터
	m_line.mx = eigen_vecs[0].x;
	m_line.my = eigen_vecs[0].y;
	m_line.sx = cntr.x;
	m_line.sy = cntr.y;
	// counter vector
//	m_line.cx = eigen_vecs[1].x;
//	m_line.cy = eigen_vecs[1].y;
	m_line.cx = -m_line.my;
	m_line.cy = m_line.mx;

	m_angle = - atan2(eigen_vecs[0].y, eigen_vecs[0].x); // orientation in radians

	Point p1(m_line.mx * -400.0 +0.5, m_line.my * -400 + 0.5);	// 0.5는 round() 대신
	Point p2(m_line.mx * 400.0 +0.5, m_line.my * 400 +0.5);

	m_p1 = cntr + p1;
	m_p2 = cntr + p2;

	return 1;
}

// PCA 직접계산
int Ransac::compute_model2(const vector<Point>& samples)
{
	// PCA(주성분분석)로 직선 백터 구한다.
	const int count = samples.size();

	double sx = 0, sy = 0;
	double sxx = 0, syy = 0, sxy = 0;

	for (int i = 0; i < count; i++)
	{
		double x = samples[i].x;
		double y = samples[i].y;

		sx += x;
		sy += y;
		sxx += x * x;
		sxy += x * y;
		syy += y * y;
	}

	//variance;
	double vxx = (sxx - sx * sx / count) / count;
	double vxy = (sxy - sx * sy / count) / count;
	double vyy = (syy - sy * sy / count) / count;

	//principal axis
	double theta = atan2(2.0 * vxy, vxx - vyy) / 2.0;

	m_line.mx = cos(theta);
	m_line.my = sin(theta);
	m_line.cx = -m_line.my;
	m_line.cy = m_line.mx;
	//center
	m_line.sx = sx / count;
	m_line.sy = sy / count;

	// 결과벡터를 두점으로 환산
	Point p1, p2;
	p1.x = round(m_line.sx - m_line.mx * 400.0);
	p1.y = round(m_line.sy - m_line.my * 400.0);
	p2.x = round(m_line.sx + m_line.mx * 400.0);
	p2.y = round(m_line.sy + m_line.my * 400.0);

	m_p1 = p1;
	m_p2 = p2;

	return 1;
}

// p1, p2를 알때
void Ransac::drawLine(Mat img, Scalar color)
{
	cv::line(img, m_p1, m_p2, color, 1);
}

void Ransac::get3PointCircle(Point2d& pt1, Point2d& pt2, Point2d& pt3, Point2d& cntr, double& r)
{
	double dy1, dy2, d, d2, yi;
	Point2d p1((pt1.x + pt2.x) / 2., (pt1.y + pt2.y) / 2.);		// 두점간 중점
	Point2d p2((pt1.x + pt3.x) / 2., (pt1.y + pt3.y) / 2.);

	dy1 = pt1.y - pt2.y;
	dy2 = pt1.y - pt3.y;
	r = 0;
	if (dy1 != 0) {
		d = (pt2.x - pt1.x) / dy1;
		yi = p1.y - d * p1.x;
		if (dy2 != 0) {
			d2 = (pt3.x - pt1.x) / dy2;
			if (d != d2) cntr.x = (yi - (p2.y - d2 * p2.x)) / (d2 - d);
			else return;
		}
		else if (pt3.x - pt1.x == 0) return;
		else cntr.x = p2.x;
	}
	else if (dy2 != 0 && pt2.x - pt1.x != 0) {
		d = (pt3.x - pt1.x) / dy2;
		yi = p2.y - d * p2.x;
		cntr.x = p1.x;
	}
	else return;
	cntr.y = d * cntr.x + yi;
	r = cv::sqrt((pt1.x - cntr.x) * (pt1.x - cntr.x) + (pt1.y - cntr.y) * (pt1.y - cntr.y));
}


bool Ransac::getBestCircle(const vector<Point>& total, const double percent, const double distance)
{
	const int total_size = total.size();
	const int max_pair = total_size * percent;
	vector<int> sample_num;
	sample_num.reserve(max_pair * 3);

	if (total.size() < 10 || max_pair < 6) return false;

	// 이전과 중복 데이터 인지 확인한다. (떨림방지)
	//if (total.size() == m_old.size())
	//{
	//	if (total == m_old) return true;
	//}

	// make random sample
	for (int i = 0; i < max_pair; i++)
	{
		int mod_num = total.size();

		int random1 = rand() % mod_num;
		int random2 = rand() % mod_num;
		int random3 = rand() % mod_num;

		if (random1 != random2 && random1 != random3 && random2 != random3)
		{
			sample_num.push_back(random1);
			sample_num.push_back(random2);
			sample_num.push_back(random3);
		}
		else i--;	// try again
	}

	// 
	int big_inliercount = 0;
	
	vector<int> sample;
	vector<int> best;
	sample.reserve(max_pair * 3);
	best.reserve(max_pair * 3);

	// RANSAC fitting
	for (int i = 0; i < sample_num.size(); i += 3)
	{
		const int p1_idx = sample_num[i];
		const int p2_idx = sample_num[i + 1];
		const int p3_idx = sample_num[i + 2];

		Point2d p[3] = { total[p1_idx], total[p2_idx], total[p3_idx] };

		//model estimation
		Point2d op;
		double r = 0;

		get3PointCircle(p[0], p[1], p[2], op, r);

		int inliercount = 0;
		sample.clear();
		// evaluation

		for (int j = 0; j < total_size; j++)
		{
			Point2d p(total[j]);
			double len = cv::sqrt((p.x - op.x) * (p.x - op.x) + (p.y - op.y) * (p.y - op.y));

			if ( fabs(r - len) < distance)
			{
				inliercount++;
				sample.push_back(j);
			}
		}

		if (big_inliercount < inliercount)
		{
			m_circle.ox = op.x;
			m_circle.oy = op.y;
			m_circle.r = r;
			big_inliercount = inliercount;
			best.swap(sample);
		}
	}

	// 최소자승법(LSM)
	Mat A(best.size(), 3, CV_64FC1);
	Mat B(best.size(), 1, CV_64FC1);
	double* ptrA = (double*)A.data;
	double* ptrB = (double*)B.data;

	for (int i = 0; i < best.size(); i++)
	{
		const int idx = best[i];
		*ptrA++ = total[idx].x;
		*ptrA++ = total[idx].y;
		*ptrA++ = 1.0;
		*ptrB++ = -(total[idx].x * total[idx].x) - (total[idx].y * total[idx].y);
	}

	// 의사역행렬을 통해 방정식의 해를 근사한다.
	cv::Mat A_pinv(3, 3, CV_64FC1);
	invert(A, A_pinv, cv::DECOMP_SVD);

	cv::Mat X = A_pinv * B;

	double fA = X.at<double>(0, 0);
	double fB = X.at<double>(1, 0);
	double fC = X.at<double>(2, 0);
	m_circle.ox = fA / -2.;
	m_circle.oy = fB / -2.;
	m_circle.r = sqrt((fA * fA) + (fB * fB) - (4 * fC)) / 2.;

	return true;
}

