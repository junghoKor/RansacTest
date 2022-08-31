// LineDetect.cpp : 이 파일에는 'main' 함수가 포함됩니다. 거기서 프로그램 실행이 시작되고 종료됩니다.
//

#include "pch.h"
#include <iostream>
#include "Ransac.h"

int g_thresole = 100;
int g_RotateAngle = 0;

Ransac ransac1;
vector<vector<Point>> contours;


class timeChecker
{
	long long before;
public:
	void start() { before = _Xtime_get_ticks(); }
	double get() { return (double)(_Xtime_get_ticks() - before) / 10000.0; }	// ms
	string getstr() { return to_string(get()); }
};


void setContours(Mat img_bin, int maxCount)
{
	Mat1b mat = img_bin;

	//	findContours(mat, contours, RETR_TREE, CHAIN_APPROX_SIMPLE);
	findContours(mat, contours, RETR_CCOMP, CHAIN_APPROX_NONE);


	// sort largest
	sort(contours.begin(), contours.end(), [](const vector<Point>& c1, const vector<Point>& c2) {
		return contourArea(c1, false) > contourArea(c2, false);
		});
	int i = 0;
	for (vector<vector<Point> >::iterator it = contours.begin(); it != contours.end();)
	{
		if (i++ > maxCount - 1) it = contours.erase(it);
		else ++it;
	}

	// 
	// Get the moments
	vector<Moments> mu(contours.size());
	vector<RotatedRect>rt(contours.size());
	vector<Point2d>center(contours.size());
	vector<Point>icenter(contours.size());
	for (size_t i = 0; i < contours.size(); i++)
	{
		mu[i] = moments(contours[i]);
		center[i].x = mu[i].m10 / mu[i].m00;
		center[i].y = mu[i].m01 / mu[i].m00;
		icenter[i].x = (int)round(center[i].x);
		icenter[i].y = (int)round(center[i].y);
		rt[i] = minAreaRect(contours[i]);
	}

	ransac1._mu = mu;
	ransac1._rt = rt;
	ransac1._center = center;
	ransac1._icenter = icenter;
}



int main()
{
	char buf[256];

	cv::Scalar colors[3];
	colors[0] = CV_RGB(255, 0, 0);
	colors[1] = CV_RGB(0, 255, 0);
	colors[2] = CV_RGB(0, 0, 255);

	namedWindow("COLOR", WINDOW_AUTOSIZE);
	namedWindow("GRAY", WINDOW_AUTOSIZE);
	namedWindow("BIN", WINDOW_AUTOSIZE);
	namedWindow("CONTROL", WINDOW_AUTOSIZE);
	resizeWindow("CONTROL", 400, 300);

	createTrackbar("tHRESOLE", "CONTROL", &g_thresole, 255, NULL);
	createTrackbar("RotateA", "CONTROL", &g_RotateAngle, 360, NULL);

	timeChecker chk;

	Mat img_org = imread("e:/temp/circlefit2.png", IMREAD_COLOR);
	Mat img_color,img_gray, img_bin;


	while (1)
	{
		img_color = img_org.clone();

		if (g_RotateAngle != 0)
		{
			Point2f center((float)(img_org.cols / 2.), (float)(img_org.rows / 2.));
			Mat r = cv::getRotationMatrix2D(center, g_RotateAngle, 1.0);
			warpAffine(img_org, img_color, r, img_org.size());
		}


		//시간체크 시작 ( 이미지 회전연산은 제외한다 )
		chk.start();
		Mat1b blur;

		cvtColor(img_color, img_gray, cv::COLOR_BGR2GRAY);
		cv::GaussianBlur(img_gray, img_gray, Size(5, 5), 0 );
		//threshold(blur, sx.img_bin, g_thresole, 255, THRESH_BINARY);
		threshold(img_gray, img_bin, g_thresole, 255, THRESH_BINARY);

		// morphology open
		Mat knel = Mat::ones(Size(7, 7), CV_8UC1);
		dilate(img_bin, img_bin, knel);
		erode(img_bin, img_bin, knel);

		//
		// START HERE ..

		Point2f cntr1;
		double r1 = 0;

		setContours(img_bin, 1);

		ransac1.getBestCircle(contours[0], 1.1, 3.0);
		cntr1.x = round(ransac1.m_circle.ox);
		cntr1.y = round(ransac1.m_circle.oy);
		r1 = ransac1.m_circle.r;

		circle(img_color, cntr1, round(r1), CV_RGB(0, 255, 100), 2);
		//putText(sx.img_color, "C1", cntr1 + Point2f(0, -90), FONT_HERSHEY_PLAIN, 1, Scalar(0, 0, 255));
		//center
		circle(img_color, cntr1, 7, CV_RGB(255,0, 0), -1);

		/// INFORMATION
		sprintf_s(buf, "C1: Radius = %.2f", r1);
		putText(img_color, buf, Point(10, 400), FONT_HERSHEY_PLAIN, 1.2, CV_RGB(255, 255, 255));
		
		imshow("GRAY", img_gray);
		imshow("BIN", img_bin);

		sprintf_s(buf, "Delay %.4fms", chk.get());
		putText(img_color, buf, Point(10, 15), FONT_HERSHEY_PLAIN, 1.2, CV_RGB(255,255,255));

		imshow("COLOR", img_color);

		int key = waitKey(10);
		if (key == 27) break;
	}
	cv::destroyAllWindows();
	return 0;
}

