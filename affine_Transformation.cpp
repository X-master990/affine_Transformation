#include<opencv2/opencv.hpp>
#include<iostream>
#include<string>
#include<fstream>
#include <iomanip>
#include <vector>
#include <math.h>
using namespace std;
using namespace cv;
#define PI 3.14159265359
bool iszero(double x, double y, double z,int l,int w ,int h) {
	if (x < 0 || x >= l)
		return true;
	if (y < 0 || y >= w)
		return true;
	if (z < 0 || z >= w)
		return true;
	return false;
}


int main(int argc, char* argv[])
{
		Mat V1 = Mat::ones(4, 1, CV_64F);
		Mat V2 = Mat::ones(4, 1, CV_64F);
		Mat V3 = Mat::ones(4, 1, CV_64F);
		Mat V4 = Mat::ones(4, 1, CV_64F);
		Mat u = Mat::ones(4, 1, CV_64F);

		fstream op(argv[1], ios::in);
		for (int i = 0; i < V1.rows-1; i++) {
			double tmp = 0;
			op >> tmp;
			V1.at<double>(i, 0) = tmp;
		}
		
		for (int i = 0; i < V2.rows-1;i++) {
			double tmp = 0;
			op >> tmp;
			V2.at<double>(i, 0) = tmp;
		}

		for (int i = 0; i < V3.rows-1; i++) {
			double tmp = 0;
			op >> tmp;
			V3.at<double>(i, 0) = tmp;
		}
		for (int i = 0; i < V4.rows-1; i++) {
			double tmp = 0;
			op >> tmp;
			V4.at<double>(i, 0) = tmp;
		}

		
		for (int i = 0; i < u.rows-1; i++) {
			double tmp = 0;
			op >> tmp;
			u.at<double>(i, 0) = tmp;
		}
		int n = 0;
		
		op >> n;
		Mat Tt = Mat::eye(4, 4, CV_64F);
		while (n--) {
			string command;
			op >> command;
			Mat T = Mat::eye(4, 4, CV_64F);

			if (command == "#T") {
				for (int i = 0; i < 3; i++) {
					double tmp=0;
					op >> tmp;
					T.at<double>(i, 3) = tmp;
				}
				
			}
			if (command == "#S"){
				double c[3];
				for (int i = 0; i < 3; i++) {
					op >> c[i];
				}
				Mat s = Mat::eye(4, 4, CV_64F);
				for (int i = 0; i < 3; i++) {
					double tmp = 0;
					op >> tmp;
					s.at<double>(i, i) = tmp;
				}
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * s * move * T;

			}
			if (command == "#Pxy") {
				T.at<double>(2, 2) = 0;
			}
			if (command == "#Pxz") {
				T.at<double>(1, 1) = 0;
			}
			if (command == "#Pyz") {
				T.at<double>(0, 0) = 0;
			}
			if (command == "#Rx") {
				double c[3];
				for (int i = 0; i < 3; i++) {
					op >> c[i];
				}
				double theata = 0;
				op >> theata;
				cv::Mat rotation_x = (cv::Mat_<double>(4, 4) << 1, 0, 0, 0,
					0, cos(theata * PI / 180), -sin(theata * PI / 180), 0,
					0, sin(theata * PI / 180), cos(theata * PI / 180), 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * rotation_x * move * T;
				
			}
			if (command == "#Ry") {
				double c[3];
				for (int i = 0; i < 3; i++) {
					op >> c[i];
				}
				double theata = 0;
				op >> theata;

				Mat rotation_y = (cv::Mat_<double>(4,4) << cos(theata * PI / 180), 0, sin(theata * PI / 180),0,
					0, 1, 0,0,
					-sin(theata * PI / 180), 0, cos(theata * PI / 180),0,
					0,0,0,1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2*rotation_y * move * T;
			}
			if (command == "#Rz") {
				
				double c[3];
				for (int i = 0; i < 3; i++) {
					op >> c[i];
				}
				double theata = 0;
				op >> theata;
				Mat rotation_z = (cv::Mat_<double>(4, 4) << cos(theata * PI / 180), -sin(theata * PI / 180), 0, 0,
					sin(theata * PI / 180), cos(theata * PI / 180), 0, 0,
					0, 0, 1,0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * rotation_z * move * T;
			}
			if (command == "#Hx") {
				double s, t;
				double c[3];
				for (int i = 0; i < 3; i++) {
					op >> c[i];
				}
				op >> s >> t;
				Mat hx= (cv::Mat_<double>(4, 4) << 1, s, t, 0,
					0,1,0, 0,
					0, 0, 1, 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * hx * move * T;
			}
			if (command == "#Hy") {
				double s, t;
				double c[3];
				for (int i = 0; i < 3; i++) {
					op >> c[i];
				}
				op >> s >> t;
				Mat hy = (cv::Mat_<double>(4, 4) << 1, 0, 0, 0,
					s, 1, t, 0,
					0, 0, 1, 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * hy * move * T;
			}
			if (command == "#Hz") {
				double s, t;
				double c[3];
				for (int i = 0; i < 3; i++) {
					op >> c[i];
				}
				op >> s >> t;
				Mat hz = (cv::Mat_<double>(4, 4) << 1, 0, 0, 0,
					0, 1, 0, 0,
					s, t, 1, 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * hz * move * T;
			}
			if (command == "#M") {
				for (int i = 0; i < T.rows ; i++) {
					for (int j = 0; j < T.cols; j++) {
						double tmp;
						op >> tmp;
						T.at<double>(i, j) = tmp;
					}
				}
			}
			Tt = T * Tt;
			
		}
		
		op.close();
		
		ofstream write(argv[2], ios::out);
		for (int i = 0; i < Tt.rows; i++) {
			for (int j = 0; j < Tt.cols; j++) {
				write << fixed << setprecision(2) << Tt.at<double>(i, j);
				if (j != Tt.cols - 1)
					write << " ";
			}
			write << endl;
		}
		//system("pause");
		Mat uV = Mat::zeros(3, 3, CV_64F);
		Mat VV = Mat::zeros(3, 3, CV_64F);
		Mat u1 = Mat::zeros(4, 1, CV_64F);
		u1=Tt * V1;
		for (int i = 0; i < u1.rows; i++) {
			for (int j = 0; j < u1.cols; j++) {
				write << fixed << setprecision(2) << u1.at<double>(i, j);
				if (i != u1.rows - 1)
					write << " ";
			}
		}
		write << endl;
		Mat u2 = Mat::zeros(4, 1, CV_64F);
		u2 = Tt * V2;
		for (int i = 0; i < u2.rows; i++) {
			for (int j = 0; j < u2.cols; j++) {
				write << fixed << setprecision(2) << u2.at<double>(i, j);
				if (i != u2.rows-1)
					write << " ";
			}
		}
		write << endl;
		Mat u3 = Mat::zeros(4, 1, CV_64F);
		u3 = Tt * V3;
		for (int i = 0; i < u3.rows; i++) {
			for (int j = 0; j < u3.cols; j++) {
				write << fixed << setprecision(2) << u3.at<double>(i, j);
				if (i != u3.rows - 1)
					write << " ";
			}
		}
		write << endl;
		Mat u4 = Mat::zeros(4, 1, CV_64F);
		u4 = Tt * V4;
		for (int i = 0; i < u4.rows; i++) {
			for (int j = 0; j < u4.cols; j++) {
				write << fixed << setprecision(2) << u4.at<double>(i, j);
				if (i != u4.rows - 1)
					write << " ";
			}
		}
		write << endl;
		Mat z = u1.clone();
		u1 = u2 - z;
		u2 = u3 - z;
		u3 = u4 - z;
		for (int i = 0; i < uV.cols; i++) {
			uV.at<double>(0, i) = u1.at<double>(i,0);
		}
		for (int i = 0; i < uV.cols; i++) {
			uV.at<double>(1, i) = u2.at<double>(i, 0);
		}
		for (int i = 0; i < uV.cols; i++) {
			uV.at<double>(2, i) = u3.at<double>(i, 0);
		}
		/*
		for (int i = 0; i < uV.cols; i++) {
			uV.at<double>(3, i) = u4.at<double>(i, 0);
		}*/
		Mat re = V1.clone();
		V1 = V2 - re;
		V2 = V3 - re;
		V3 = V4 - re;
		for (int i = 0; i < VV.cols; i++) {
			VV.at<double>(0, i) = V1.at<double>(i, 0);
		}
		for (int i = 0; i < VV.cols; i++) {
			VV.at<double>(1, i) = V2.at<double>(i, 0);
		}
		for (int i = 0; i < VV.cols; i++) {
			VV.at<double>(2, i) = V3.at<double>(i, 0);
		}
		/*
		for (int i = 0; i < VV.cols; i++) {
			VV.at<double>(3, i) = V4.at<double>(i, 0);
		}*/


		double B = determinant(uV);
		double A = determinant(VV);
		
		
		write << fixed << setprecision(2) << abs(B / A) << " " << determinant(Tt) << endl;
		
		float r = abs(B / A);
		float dt = determinant(Tt);
		//cout << (r == dt) << endl;
		
		
		if (int(r) == 0 && int(dt) == 0) {
			write << "zeros" << endl;
		}
		else if (r == -dt) {
			write << "r==-det(T)" << endl;
		}
		else if (r == dt) {
			write << "r==det(T)" << endl;
		}
		else {
			write << "others" << endl;
		}

		Mat Tt_inverse = Mat::zeros(4, 4, CV_64F);
		invert(Tt, Tt_inverse, 0);
		Mat v = Mat::zeros(4, 4, CV_64F);
		v = Tt_inverse * u;
		if (determinant(Tt) == 0) {
			write << "NaN";
		}
		else
			for (int i = 0; i < v.rows; i++) {
				for (int j = 0; j < v.cols; j++) {
					write << v.at<double>(i, j);
					if (i != v.rows - 1)
						write << " ";
				}
		}
		write << endl;
		write.close();
		
		ifstream open(argv[3], ios::in);
		int l, w, h;
		open >> l >> w >> h;
		
		vector<Mat> data;
		
		while (h--) {
			Mat tmp = Mat(l, w, CV_64F);
			for (int i = 0; i < l; i++) {
				for (int j = 0; j < w; j++) {
					open >> tmp.at<double>(i, j);
				}
			}
			data.push_back(tmp);
		}
		int k = 0;
		open >> k;
		Mat tr = Mat::eye(4, 4, CV_64F);
		while (k--) {
			string command;
			open >> command;
			Mat T = Mat::eye(4, 4, CV_64F);

			if (command == "#T") {
				for (int i = 0; i < 3; i++) {
					double tmp = 0;
					open >> tmp;
					T.at<double>(i, 3) = tmp;
				}

			}
			if (command == "#S") {
				double c[3];
				for (int i = 0; i < 3; i++) {
					open >> c[i];
				}
				Mat s = Mat::eye(4, 4, CV_64F);
				for (int i = 0; i < 3; i++) {
					double tmp = 0;
					open >> tmp;
					s.at<double>(i, i) = tmp;
				}
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * s * move * T;

			}
			if (command == "#Pxy") {
				T.at<double>(2, 2) = 0;
			}
			if (command == "#Pxz") {
				T.at<double>(1, 1) = 0;
			}
			if (command == "#Pyz") {
				T.at<double>(0, 0) = 0;
			}
			if (command == "#Rx") {
				double c[3];
				for (int i = 0; i < 3; i++) {
					open >> c[i];
				}
				double theata = 0;
				open >> theata;
				cv::Mat rotation_x = (cv::Mat_<double>(4, 4) << 1, 0, 0, 0,
					0, cos(theata * PI / 180), -sin(theata * PI / 180), 0,
					0, sin(theata * PI / 180), cos(theata * PI / 180), 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * rotation_x * move * T;

			}
			if (command == "#Ry") {
				double c[3];
				for (int i = 0; i < 3; i++) {
					open >> c[i];
				}
				double theata = 0;
				open >> theata;

				Mat rotation_y = (cv::Mat_<double>(4, 4) << cos(theata * PI / 180), 0, sin(theata * PI / 180), 0,
					0, 1, 0, 0,
					-sin(theata * PI / 180), 0, cos(theata * PI / 180), 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * rotation_y * move * T;
			}
			if (command == "#Rz") {

				double c[3];
				for (int i = 0; i < 3; i++) {
					open >> c[i];
				}
				double theata = 0;
				open >> theata;
				Mat rotation_z = (cv::Mat_<double>(4, 4) << cos(theata * PI / 180), -sin(theata * PI / 180), 0, 0,
					sin(theata * PI / 180), cos(theata * PI / 180), 0, 0,
					0, 0, 1, 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * rotation_z * move * T;
			}
			if (command == "#Hx") {
				double s, t;
				double c[3];
				for (int i = 0; i < 3; i++) {
					open >> c[i];
				}
				open >> s >> t;
				Mat hx = (cv::Mat_<double>(4, 4) << 1, s, t, 0,
					0, 1, 0, 0,
					0, 0, 1, 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * hx * move * T;
			}
			if (command == "#Hy") {
				double s, t;
				double c[3];
				for (int i = 0; i < 3; i++) {
					open >> c[i];
				}
				open >> s >> t;
				Mat hy = (cv::Mat_<double>(4, 4) << 1, 0, 0, 0,
					s, 1, t, 0,
					0, 0, 1, 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * hy * move * T;
			}
			if (command == "#Hz") {
				double s, t;
				double c[3];
				for (int i = 0; i < 3; i++) {
					open >> c[i];
				}
				open >> s >> t;
				Mat hz = (cv::Mat_<double>(4, 4) << 1, 0, 0, 0,
					0, 1, 0, 0,
					s, t, 1, 0,
					0, 0, 0, 1);
				Mat move = (cv::Mat_<double>(4, 4) << 1, 0, 0, -c[0],
					0, 1, 0, -c[1],
					0, 0, 1, -c[2],
					0, 0, 0, 1);
				Mat move2 = (cv::Mat_<double>(4, 4) << 1, 0, 0, c[0],
					0, 1, 0, c[1],
					0, 0, 1, c[2],
					0, 0, 0, 1);
				T = move2 * hz * move * T;
			}
			if (command == "#M") {
				for (int i = 0; i < T.rows; i++) {
					for (int j = 0; j < T.cols; j++) {
						double tmp;
						open >> tmp;
						T.at<double>(i, j) = tmp;
					}
				}
			}
			tr = T * tr;
			
		}
		Mat tr_inverse = Mat::zeros(4, 4, CV_64F);
		invert(tr, tr_inverse, 0);
	
		ofstream add(argv[4], ios::out);

		h = data.size();
		

		for (int i = 0; i < data.size(); i++) {
			for (int j = 0; j < data[i].rows; j++) {
				for (int k = 0; k < data[i].cols; k++) {
					Mat p = (Mat_<double>(4, 1) <<
						j,
						k,
						i,
						1);
					
					Mat A = tr_inverse * p;

					int P[3] = {floor(A.at<double>(0, 0)),
						floor(A.at<double>(1, 0)),
						floor(A.at<double>(2, 0))};
					//cout << P[0] << " " << P[1] << " " << P[2] << endl;
					int arr[2][2][2];
					for (int a = 0; a < 2; a++) {
						for (int b = 0; b < 2; b++) {
							for (int c = 0; c < 2; c++) {
								if (a+P[0] < 0 || a+P[0] >= l || b+P[1] < 0 || b+P[1] >= w || c+P[2] < 0 || c+P[2] >= h) {
									arr[a][b][c] = 0;
								}
								else {
									arr[a][b][c] = data[c+P[2]].at<double>(a+P[0], b+P[1]);
								}
							}
						}
					}
					
					/*
					double t = a.at<double>(0, 0);
					double t1 = t;
					t = floor(t);
					t1 = ceil(t1);
					double x0 = t;
					double x1 = t1;

					t = a.at<double>(1, 0);
					t1 = t;
					t = floor(t);
					t1 = ceil(t1);
					double y0 = t;
					double y1 = t1;
					t = a.at<double>(2, 0);
					t1 = t;
					t = floor(t);
					t1 = ceil(t1);
					double z0 = t;
					double z1 = t1;

					double xd;
					if ((x1 - x0) <= 0)
						xd = 0;
					else
						xd = (a.at<double>(0, 0) - x0) / (x1 - x0);
					double yd;
					if ((y1 - y0) <= 0)
						yd = 0;
					else
						yd = (a.at<double>(1, 0) - y0) / (y1 - y0);
					double zd;
					if ((z1 - z0) <= 0)
						zd = 0;
					else
						zd = (a.at<double>(2, 0) - z0) / (z1 - z0);
					//cout << x1 << "and" << y1 << "and" << z1 << endl;
					//cout << x0 << "and" << y0 << "and" << z0 << endl;
					//cout << xd << "and" << yd << "and" << zd << endl;

					double c_00;
					if (z0 < 0 || x0 < 0 || x0 >= l || z0 >= h || (y0 >= w && y1 >= w))
						c_00 = 0;
					else if (y0 < 0 || y0 >= w)
						c_00 = data[z0].at<double>(x0, y1) * xd;
					else if (y1 < 0 || y1 >= w)
						c_00 = data[z0].at<double>(x0, y0) * (1 - xd);
					else {
						//cout << x0 << "and" << y0 << "and" << z0 << endl;
						c_00 = data[z0].at<double>(x0, y0) * (1 - xd) + data[z0].at<double>(x0, y1) * xd;
					}

					double c_01;

					if (z1 < 0 || x0 < 0 || x0 >= l || z1 >= h || (y1 >= w && y0 >= w))
						c_01 = 0;
					else if (y0 < 0 || y0 >= w)
						c_01 = data[z1].at<double>(x0, y1) * xd;
					else if (y1 < 0 || y1 >= w)
						c_01 = data[z1].at<double>(x0, y0) * (1 - xd);
					else {
						//cout << x0 << "and" << y0 << "and" << z0 << endl;
						c_01 = data[z1].at<double>(x0, y0) * (1 - xd) + data[z1].at<double>(x0, y1) * xd;
					}
					//cout << "c01:" << c_01 << endl;

					double c_10;
					if (z0 < 0 || x1 < 0 || x1 >= l || z0 >= h || (y1 >= w && y0 >= w))
						c_10 = 0;
					else if (y0 < 0 || y0 >= w)
						c_10 = data[z0].at<double>(x1, y1) * xd;
					else if (y1 < 0 || y1 >= w)
						c_10 = data[z0].at<double>(x1, y0) * (1 - xd);
					else {
						//cout << x0 << "and" << y0 << "and" << z0 << endl;
						c_10 = data[z0].at<double>(x1, y0) * (1 - xd) + data[z0].at<double>(x1, y1) * xd;
					}

					double c_11;
					if (z1 < 0 || x1 < 0 || x1 >= l || z1 >= h || (y0 >= w && y1 >= w))
						c_11 = 0;
					else if (y0 < 0 || y0 >= w)
						c_11 = data[z1].at<double>(x1, y1) * xd;
					else if (y1 < 0 || y1 >= w)
						c_11 = data[z1].at<double>(x1, y0) * (1 - xd);
					else {
						//cout << x0 << "and" << y0 << "and" << z0 << endl;
						c_11 = data[z1].at<double>(x1, y0) * (1 - xd) + data[z1].at<double>(x1, y1) * xd;
					}
					//cout << c_00 << "and" << c_01 << "and" << c_10 << "and" << c_11 << endl;

*/
					double X[4];
					double Y[2];
					int C;
					X[0] = arr[1][0][1] * (A.at<double>(0, 0) - P[0]) + arr[0][0][1] * (P[0] + 1 - A.at<double>(0, 0));
					X[1] = arr[1][1][1] * (A.at<double>(0, 0) - P[0]) + arr[0][1][1] * (P[0] + 1 - A.at<double>(0, 0));
					X[2] = arr[1][0][0] * (A.at<double>(0, 0) - P[0]) + arr[0][0][0] * (P[0] + 1 - A.at<double>(0, 0));
					X[3] = arr[1][1][0] * (A.at<double>(0, 0) - P[0]) + arr[0][1][0] * (P[0] + 1 - A.at<double>(0, 0));
					Y[0] = X[0] * (P[1] + 1 - A.at<double>(1, 0)) + X[1] * (A.at<double>(1, 0) - P[1]);
					Y[1] = X[2] * (P[1] + 1 - A.at<double>(1, 0)) + X[3] * (A.at<double>(1, 0) - P[1]);
					C = Y[0] * (A.at<double>(2, 0) - P[2]) + Y[1] * (P[2] + 1 - A.at<double>(2, 0));
					add << C;
					if (k != data[i].cols - 1)
						add << " ";
				}
				add << endl;
			}
			
		}
		

		
}