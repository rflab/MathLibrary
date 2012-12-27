//unit test
//rflab.

#include <iostream>
#include <iomanip>

//rock field lib
#include "math/index.h"
#include "scopetime.h"

//DirectX
#pragma comment(lib, "d3dx9.lib")
#include <d3dx9.h>

//boost::ublas
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


int main()
{
	using namespace rf::math;
	using namespace std;

	//basic math
	#if 1
		{
			cout << "\n>>math\n" << endl;
			cout << function::RandUnit()        << endl;
			cout << function::RandUnit()        << endl;
			cout << function::RandUnit()        << endl;
			cout << function::Sqrt(64)          << endl;
			cout << function::RSqrt(64)         << endl;
			cout << function::Sin(3.14f)        << endl;
			cout << function::Cos(3.14f)        << endl;
			cout << function::Tan(3.14f)        << endl;
			cout << function::ASin(1.0f)        << endl;
			cout << function::ACos(1.0f)        << endl;
			cout << function::ATan2(1.0f, 1.0f) << endl;
		}
	#endif

	//quaternion
	#if 1
		//生成
		#if 1
			{
				cout << "\n>>transform\n" << endl;

				t_quaternion q1(1, 2, 3, 4);
				cout << (string)q1 << endl;
				
				t_quaternion q2;
				
				function::MakeQuaternionIdentity(&q2);
				cout << (string)q2 << endl;

				function::MakeQuaternionNormalize(q1, &q2);
				cout << (string)q2 << endl;
				cout << q2.Length() << endl;
				
				function::MakeQuaternionInverse(q1, &q2);
				cout << (string)q2 << endl;
				cout << (string)(q1*q2) << endl;

				function::MakeQuaternionConjugate(q1, &q2);
				cout << (string)q2 << endl;
				
				function::MakeQuaternionRotationX(definition::RF_PI_OVER2, &q2);
				cout << (string)q2 << endl;
				function::MakeQuaternionRotationY(definition::RF_PI_OVER2, &q2);
				cout << (string)q2 << endl;
				function::MakeQuaternionRotationY(definition::RF_PI_OVER2, &q2);
				cout << (string)q2 << endl;
			}
		#endif

		//球面線形補間
		#if 1
			{
				cout << ">>slerp" << endl;

				t_quaternion q1(1, 2, -3, -4);
				t_quaternion q2(1, -1, -3, 3);
				t_quaternion q3;
				q1.Normalize();
				q2.Normalize();
				
				t_quaternion c;

				for (float f = 0; f<=1.0f; f+=0.1f)
				{
					function::MakeQuaternionSlerp(q1, q2, f, &q3);
					cout << "amount("<< f << ")" << (string)c << endl;
				}
			}
		#endif
		
		//変換
		#if 1
			{
				//Y軸で90度回転するクォータニオンから変換
				cout << "\n>>quaternion\n" << endl;

				t_quaternion q;
				t_euler_angles e;
				t_matrix4 m;

				function::MakeQuaternionRotationY(definition::RF_PI_OVER2, &q);
				function::MakeMatrixFromQuaternion(q, &m);
				function::MakeEulerAnglesFromQuaternion(q, &e);

				cout << (string)q << endl;
				cout << (string)e << endl;
				cout << (string)m << endl;
			}

			{
				//Y軸で90度回転する行列から変換
				cout << "\n>>rotation matrix\n" << endl;

				t_quaternion q;
				t_euler_angles e;
				t_matrix4 m;
			
				function::MakeMatrixRotationY(definition::RF_PI_OVER2, &m);
				function::MakeQuaternionFromMatrix(m, &q);
				function::MakeEulerAnglesFromMatrix(m, &e);

				cout << (string)q << endl;
				cout << (string)e << endl;
				cout << (string)m << endl;
			}

			{
				//Y軸で90度回転するオイラー角から変換
				cout << "\n>>euler angles\n" << endl;

				t_euler_angles e(3.141592f/2.0f, 0, 0);
				t_quaternion q;
				t_matrix4 m;

				function::MakeQuaternionFromEulerAngles(e, &q);
				function::MakeMatrixFromEulerAngles(e, &m);

				cout << (string)q << endl;
				cout << (string)e << endl;
				cout << (string)m << endl;
			}
		#endif
			
		//回転
		#if 1
			{
				//クォータニオンでベクトルを回転
				cout << "\n>>quaternion rotate\n" << endl;
				t_quaternion qx;
				t_quaternion qy;
				t_quaternion qz;
				function::MakeQuaternionFromAxisAngle(t_vector3(1, 0, 0), 3.141592f/2.0f, &qx);
				function::MakeQuaternionFromAxisAngle(t_vector3(0, 1, 0), 3.141592f/2.0f, &qy);
				function::MakeQuaternionFromAxisAngle(t_vector3(0, 0, 1), 3.141592f/2.0f, &qz);

				t_vector3 a(1, 1, 1);
				function::Rotate(qx, &a);
				cout << "x rot "<< (string)a << endl;

				a = t_vector3(1, 1, 1);
				function::Rotate(qy, &a);
				cout << "y rot "<< (string)a << endl;

				a = t_vector3(1, 1, 1);
				function::Rotate(qz, &a);
				cout << "z rot "<< (string)a << endl;
			}
		
			//クォータニオンを行列に変換してベクトルを回転
			{
				cout << "\n>>quaternion to matrix and rotate\n" << endl;
				t_quaternion qx;
				t_quaternion qy;
				t_quaternion qz;
				function::MakeQuaternionFromAxisAngle(t_vector3(1, 0, 0), 3.141592f/2.0f, &qx);
				function::MakeQuaternionFromAxisAngle(t_vector3(0, 1, 0), 3.141592f/2.0f, &qy);
				function::MakeQuaternionFromAxisAngle(t_vector3(0, 0, 1), 3.141592f/2.0f, &qz);
				
				//クォータニオンを行列に変換
				t_matrix4 mx;
				t_matrix4 my;
				t_matrix4 mz;
				function::MakeMatrixFromQuaternion(qx, &mx);
				function::MakeMatrixFromQuaternion(qy, &my);
				function::MakeMatrixFromQuaternion(qz, &mz);

				t_vector3 a(1, 1, 1);
				function::TransformCoordinate(mx, &a);
				cout << "x rot "<< (string)a << endl;

				a = t_vector3(1, 1, 1);
				function::TransformCoordinate(my, &a);
				cout << "y rot "<< (string)a << endl;

				a = t_vector3(1, 1, 1);
				function::TransformCoordinate(mz, &a);
				cout << "z rot "<< (string)a << endl;
			}
		#endif
	#endif

	//casting
	#if 1
		{
			cout << "\n>>casting\n" << endl;
			D3DXVECTOR4 d3dVector(
				1,2,3,4);

			D3DXMATRIX d3dMatrix(
				 1, 2, 3, 4,
				 5, 6, 7, 8,
				 9,10,11,12,
				13,14,15,16);

			float fVector[] = {
				0,0,0,0};

			float fMatrix[] = {
				0,0,0,0,
				0,0,0,0,
				0,0,0,0,
				0,0,0,0};

			{
				//d3d copy constructor
				t_vector4 v(d3dVector);
				t_matrix4 m(d3dMatrix);
				cout << (string)v << endl;
				cout << (string)m << endl;
			}

			{
				//float array copy constructor
				t_vector4 v(fVector);
				t_matrix4 m(fMatrix);
				cout << (string)v << endl;
				cout << (string)m << endl;
			}

			{
				//d3d copy
				t_vector4 v;
				t_matrix4 m;
				v = d3dVector;
				m = d3dMatrix;
				cout << (string)v << endl;
				cout << (string)m << endl;
			}

			{
				//float array copy
				t_vector4 v;
				t_matrix4 m;
				v = fVector;
				m = fMatrix;
				cout << (string)v << endl;
				cout << (string)m << endl;
			}
		}
	#endif


	//ベクトル速度比較
	#if 1
		#if 1
			//rf vector
			{
				cout << "my vector4" << endl;
				float _b11, _b12, _b13, _b14;
				float _a11, _a12, _a13, _a14;

				cout<< "matrix a input =" << endl;
				cin >> _a11 >> _a12 >> _a13 >> _a14;

				cout<< "matrix b input =" << endl;
				cin >> _b11 >> _b12 >> _b13 >> _b14;

				t_vector4 a;
				a(0) = _a11;
				a(1) = _a12;
				a(2) = _a13;
				a(3) = _a14;
							
									
				cout << (string)a << endl;

				int loop;
				for(;;)
				{
					t_vector4 b;
					b(0) = _b11;
					b(1) = _b12;
					b(2) = _b13;
					b(3) = _b14;

					cin >> loop;
					if (cin.fail())
						break;

					cout<< "loop";
					rf::CClock clk;
					for(int i=0;i<loop;i++)
					{
						a = a+b+b+b+b+b+b;
					}
					cout << (string)clk << endl;

					cout << (string)(b) << endl;
					cout << (string)(a)<< endl;
				}
				cin.clear();
				cin.seekg(0);
				
			}
		#elif 1
			//d3d vector
			{
				cout << "d3d" << endl;
				float _b11, _b12, _b13, _b14;
				float _a11, _a12, _a13, _a14;

				cout<< "matrix a input =" << endl;
				cin >> _a11 >> _a12 >> _a13 >> _a14;

				cout<< "matrix b input =" << endl;
				cin >> _b11 >> _b12 >> _b13 >> _b14;


				D3DXVECTOR4 a(_a11, _a12, _a13, _a14);
							
						
				cout << a.x << a.y << a.z << a.w << endl;

				int loop;
				for(;;)
				{
					D3DXVECTOR4 b(_b11, _b12, _b13, _b14);

					cout<< "loop";
					cin >> loop;
					if (cin.fail())
						break;

					rf::CClock clk;
					for(int i=0;i<loop;i++)
					{
						a = a+b+b+b+b+b+b;
					}
					cout << (string)clk << endl;

					cout << a.x << a.y << a.z << a.w << endl;
					cout << b.x << b.y << b.z << b.w << endl;
				}
			}
		#endif
	#endif

	//行列速度比較
	#if 1
		#if 1
			{
				cout << "my matrix" << endl;
				float _b11, _b12, _b13, _b14,
					  _b21, _b22, _b23, _b24,
					  _b31, _b32, _b33, _b34,
					  _b41, _b42, _b43, _b44;

				float _a11, _a12, _a13, _a14,
					  _a21, _a22, _a23, _a24,
					  _a31, _a32, _a33, _a34,
					  _a41, _a42, _a43, _a44;

				cout<< "matrix a input =" << endl;
					cin >> _a11 >> _a12 >> _a13 >> _a14 >>
					_a21 >> _a22 >> _a23 >> _a24 >>
					_a31 >> _a32 >> _a33 >> _a34 >>
					_a41 >> _a42 >> _a43 >> _a44;

				cout<< "matrix b input =" << endl;
				cin >> _b11 >> _b12 >> _b13 >> _b14 >>
					_b21 >> _b22 >> _b23 >> _b24 >>
					_b31 >> _b32 >> _b33 >> _b34 >>
					_b41 >> _b42 >> _b43 >> _b44;

				t_matrix4 a(
					_a11, _a12, _a13, _a14,
					_a21, _a22, _a23, _a24,
					_a31, _a32, _a33, _a34,
					_a41, _a42, _a43, _a44);
						
				cout << (string)a << endl;

				int loop;
				for(;;)
				{
					t_matrix4 b(
						_b11, _b12, _b13, _b14,
						_b21, _b22, _b23, _b24,
						_b31, _b32, _b33, _b34,
						_b41, _b42, _b43, _b44);

					t_vector3 c(_a11, _a12, _a13);
					t_vector3 d(_b11, _b12, _b13);

					cout << "loop n = ";
					cin  >> loop;
					if (cin.fail())
						break;

					rf::CClock clk;
					for(int i=0;i<loop;i++)
					{
						//b *= a;
						//b = d+c+c+c+c+c+c;
						b = b*a;
					}
					cout << (string)clk << endl;

					cout << (string)b << endl;
				}

				cin.clear();
				cin.seekg(0);
				return 0;
			}
		#elif 1
			{
				cout << "d3d" << endl;
				float _b11, _b12, _b13, _b14,
					  _b21, _b22, _b23, _b24,
					  _b31, _b32, _b33, _b34,
					  _b41, _b42, _b43, _b44;

				float _a11, _a12, _a13, _a14,
					  _a21, _a22, _a23, _a24,
					  _a31, _a32, _a33, _a34,
					  _a41, _a42, _a43, _a44;

				cout<< "matrix a input =" << endl;
				cin >> _a11 >> _a12 >> _a13 >> _a14 >>
					_a21 >> _a22 >> _a23 >> _a24 >>
					_a31 >> _a32 >> _a33 >> _a34 >>
					_a41 >> _a42 >> _a43 >> _a44;

				cout<< "matrix b input =" << endl;
				cin >> _b11 >> _b12 >> _b13 >> _b14 >>
					_b21 >> _b22 >> _b23 >> _b24 >>
					_b31 >> _b32 >> _b33 >> _b34 >>
					_b41 >> _b42 >> _b43 >> _b44;


				D3DXMATRIX a(
					_a11, _a12, _a13, _a14,
					_a21, _a22, _a23, _a24,
					_a31, _a32, _a33, _a34,
					_a41, _a42, _a43, _a44);
							
						
				cout << a._11 << a._12 << a._13 << a._14 << endl;
				cout << a._21 << a._22 << a._23 << a._24 << endl;
				cout << a._31 << a._32 << a._33 << a._34 << endl;
				cout << a._41 << a._42 << a._43 << a._44 << endl;

				int loop;
				for(;;)
				{
					D3DXMATRIX b(
						_b11, _b12, _b13, _b14,
						_b21, _b22, _b23, _b24,
						_b31, _b32, _b33, _b34,
						_b41, _b42, _b43, _b44);

					D3DXVECTOR3 c(_a11, _a12, _a13);
					D3DXVECTOR3 d(_b11, _b12, _b13);

					cout << "loop n = ";
					cin  >> loop;
					if (cin.fail())
						break;

					rf::CClock clk;
					for(int i=0;i<loop;i++)
					{
						b = b * a;
					}
					cout << (string)clk << endl;

					cout << 	b._11 << b._12 << b._13 << b._14 << endl;
					cout << 	b._21 << b._22 << b._23 << b._24 << endl;
					cout << 	b._31 << b._32 << b._33 << b._34 << endl;
					cout << 	b._41 << b._42 << b._43 << b._44 << endl;
				}

				cin.clear();
				cin.seekg(0);
				return 0;
			}
		#elif 1
			{
				using namespace boost::numeric::ublas;

				float aArray[4][4];
				float bArray[4][4];

				cout<< "matrix a input =" << endl;
				cin >> aArray[1][1] >> aArray[1][2] >> aArray[1][3] >> aArray[1][4] >>
					aArray[2][1] >> aArray[2][2] >> aArray[2][3] >> aArray[2][4] >>
					aArray[3][1] >> aArray[3][2] >> aArray[3][3] >> aArray[3][4] >>
					aArray[4][1] >> aArray[4][2] >> aArray[4][3] >> aArray[4][4];
				
				cout<< "matrix b input =" << endl;
				cin >> bArray[1][1] >> bArray[1][2] >> bArray[1][3] >> bArray[1][4] >>
					bArray[2][1] >> bArray[2][2] >> bArray[2][3] >> bArray[2][4] >>
					bArray[3][1] >> bArray[3][2] >> bArray[3][3] >> bArray[3][4] >>
					bArray[4][1] >> bArray[4][2] >> bArray[4][3] >> bArray[4][4];

				matrix<float> a( 4, 4 );
				for( int i=0; i!=4; ++i )
				for( int j=0; j!=4; ++j )
				{
					a( i, j ) = aArray[i][j];
				}
				
				cout << a << endl;


				int loop;
				for(;;)
				{	
					matrix<float> b( 4, 4 );
					for( int i=0; i!=4; ++i )
					for( int j=0; j!=4; ++j )
					{
						b( i, j ) = bArray[i][j];
					}

					cout << "loop n = ";
					cin  >> loop;
					if (cin.fail())
						break;

					rf::CClock clk;
					for(int i=0;i<loop;i++)
					{
						b = prod(b, a);
					}
					cout << (string)clk << endl;

					cout << b << endl;
				}

				cin.clear();
				cin.seekg(0);
				return 0;
			}
		#endif
	#endif
}
