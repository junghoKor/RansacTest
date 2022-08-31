#pragma once


struct sLine {
	double mx, my;
	double sx, sy;
	double cx, cy;
};

struct sCircle {
	double ox, oy;
	double r;
};




class Ransac
{
private:
	vector<Point> m_old;
	double m_a;
	double m_b;
	double m_c;

public:
	float m_score;
	vector<Point> m_best;
	sLine m_line;
	sCircle m_circle;
	Point m_p1;
	Point m_p2;
	double m_angle;
	Mat m_best3;

	vector<Moments> _mu;
	vector<RotatedRect> _rt;
	vector<Point2d> _center;
	vector<Point> _icenter;

	Ransac();
	void Clear();

	bool getBestModel(const vector<Point>& total, const double percent, const double distance);
	void drawLine(Mat img, Scalar color);
	bool getBestModel2(const vector<Point>& total, const double percent, const double distance);
	void drawLine2(Mat img, Scalar color);
	bool getBestModel3(const vector<Point>& total, const double percent, const double distance);
	void drawCurve(Mat img, Scalar color);


	int compute_model(const vector<Point>& samples);	// PCA(주성분분석)
	int compute_model2(const vector<Point>& samples);

	void get3PointCircle(Point2d& pt1, Point2d& pt2, Point2d& pt3, Point2d& cntr, double& r);
	bool getBestCircle(const vector<Point>& total, const double percent, const double distance);
};

