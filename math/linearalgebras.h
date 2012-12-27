//2010/04/18
//rflab.

#ifndef _RF_LINEAR_ALGEBLAS_FUNCTIONS_
#define _RFLINEAR_ALGEBLAS_FUNCTIONS_

#include "math.h"
#include "expressiontemplate.h"
#include "vectortemplate.h"
#include "matrixtemplate.h"
#include "euleranglestemplate.h"
#include "quaterniontemplate.h"
#include "expressionoperators.h"

//debug
#undef LOCALPRINTF
//#include <stdio.h>
//#define LOCALPRINTF(...) printf(__VA_ARGS__)
#define LOCALPRINTF(...)


///rock field!
namespace rf
{
	///数学
	namespace math
	{
		//-------------------------------------------------
		//明示的な実体化
		typedef TVector<float, 2> t_vector2;
		typedef TVector<float, 3> t_vector3;
		typedef TVector<float, 4> t_vector4;

		typedef TMatrix<float, 2, 2> t_matrix2;
		typedef TMatrix<float, 3, 3> t_matrix3;
		typedef TMatrix<float, 4, 4> t_matrix4;

		typedef TEulerAngles<float> t_euler_angles;

		typedef TQuaternion<float> t_quaternion;


		//-------------------------------------------------
		//単位行列等をグローバルstaticにおいておくと便利
		static const t_vector2 g_vector2_zero(0.0f, 0.0f);
		static const t_vector3 g_vector3_zero(0.0f, 0.0f, 0.0f);
		static const t_vector4 g_vector4_zero(0.0f, 0.0f, 0.0f, 0.0f);

		static const t_matrix2 g_matrix2_zero(
			0.0f, 0.0f,
			0.0f, 0.0f);
		static const t_matrix3 g_matrix3_zero(
			0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f);
		static const t_matrix4 g_matrix4_zero(
			0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 0.0f);

		static const t_matrix2 g_matrix2_identity(
			1.0f, 0.0f,
			0.0f, 1.0f);
		static const t_matrix3 g_matrix3_identity(
			1.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 1.0f);
		static const t_matrix4 g_matrix4_identity(
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f);
		
		static const t_quaternion g_quaternion_identity(
			0.0f, 0.0f, 0.0f, 1.0f);

		static const t_euler_angles g_euler_angles_identity(
			0.0f, 0.0f, 0.0f);


		//-------------------------------------------------for Vector
		///グローバル関数
		namespace function
		{

			///内積               <br>
			///<br>
			///a・b = ∥a∥∥b∥cosθ <br>
			///<br>
			template<typename T, int Dim>
			T Dot(const TVector<T, Dim>& lhs_, const TVector<T, Dim>& rhs_)
			{
				T tmp = m_elem[0] * rhs_.m_elem[0];
				for (int i = 1; i < Dim; ++i)
					tmp += (m_elem[i] * rhs_.m_elem[i]);
				return tmp;
			}

            ///3次元外積                                   <br>
            ///<br>
            ///二つのベクトルに直交したベクトル(トルク)を得る      <br>
            ///得られたベクトルの大きさは以下                      <br>
            ///<br>
            ///∥a×b∥ = ∥a∥∥b∥sinθ                      <br>
            ///<br>
            ///ちなみに、2次元外積は本来定義されていない           <br>
            ///<br>
            ///内積はn次数ベクトル空間で一般に成立します。     <br>
            ///（n：任意の自然数）                             <br>
            ///<br>
            ///外積はuとvの対(u,v)からEへの双写像で、          <br>
            ///u×v = -v×u                                <br>
            ///(u×v・w〉=〈u・v×w〉                      <br>
            ///を満たすものと定義します。                      <br>
            ///外積は空間の次元が1, 3, 7の時だけ定義されます。 <br>
            ///<br>
            ///上の議論はn次数ベクトル空間での話ですが、       <br>
            ///通常のユークリッド空間は数ベクトル空間ですので、<br>
            ///同じ議論で大丈夫です。                          <br>
            ///<br>
            ///とのこと。                                          <br>
            ///<br>
			template<typename T>
			inline TVector<T, 3> Cross(
				const TVector<T, 3>& lhs_,
				const TVector<T, 3>& rhs_)
			{
				return TVector<T, 3>(
					obj[0] = lhs_[1]*rhs_[2] - lhs_[2]*rhs_[1],
					obj[1] = lhs_[2]*rhs_[0] - lhs_[0]*rhs_[2],
					obj[2] = lhs_[0]*rhs_[1] - lhs_[1]*rhs_[0]);
			}
			
			///スカラ三重積                           <br>
			///<br>
			///u・v×w = ∥u∥∥v×w∥cosθ               <br>
			///<br>
			///・符号付き体積を得る                           <br>
			///・3つのベクトルの回転方向を判定できる          <br>
			///・行列式determinant([u;v;w])と等しい           <br>
			///・以下の項等式が成り立つ                       <br>
			///u・v×w = v・w×u = w・u×v                <br>
			///= -(w・v×u) = -(v・u×w) = -(u・w×v) <br>
			///<br>
			template<typename T, int Dim>
			inline T ScalarTriple(
					const TVector<T, Dim>& u_,
					const TVector<T, Dim>& v_,
					const TVector<T, Dim>& w_)
			{
				return Dot(w_, Cross(u_, v_));
			}

			///ベクトル三重積                    <br>
			///<br>
			///u×(v×w) = (w・u)v - (u・v)w         <br>
			///<br>
			///・vとwの張る平面に乗るベクトルが得られる  <br>
			///・以下の項等式が成り立つ                  <br>
			///u×(v×w) + v×(w×u) + w×(u×v) = 0 <br>
			///・何に使うのか分かりません                <br>
			///<br>
			template<typename T, int Dim>
			inline T VectorTriple(
					const TVector<T, Dim>& u_,
					const TVector<T, Dim>& v_,
					const TVector<T, Dim>& w_)
			{
				return Cross(u_, Cross(u_, v_));
			}

			///同次変換行列の回転成分のみをベクトルに適用
			template<typename T>
			inline void TransformCoordinate(
					const TMatrix<T, 4, 4>& m_,
					TVector<T, 3>* v_)
			{
				TVector<T, 3> org(*v_);

				for (int i = 0; i < 3; ++i)
				{
					(*v_)(i) = org(0) * m_(0, i);
					for (int r = 1; r < 3; ++r)
						(*v_)(i) += org(r) * m_(r, i);
				}
			}

			///クォータニオンで回転                                            <br>
			///<br>
			///q*p*(q^-1)を展開し、任意軸周りの回転行列と比較すると証明できるそう。<br>
			///幾何学的には同じ方向に角度を二分割して回転しつつ                    <br>
			///実軸要素の項等を保つような感じ？                                    <br>
			///<br>
			template<template<typename>class Quaternion, typename T>
			inline void Rotate(const Quaternion<T>& q_, TVector<T, 3>* out_)
			{
				Quaternion<T> pos;
				Quaternion<T> iq(q_);
				
				//回転したいベクトルを位置を表すクォータニオンに変換
				function::MakeQuaternionPure(*out_, &pos);

				//回転クォータニオンの逆
				iq.Inverse();

				//回転処理
				Quaternion<T> tmp(q_ * pos * iq);

				//回転後のクォータニオンをベクトルに変換
				(*out_)[0] = tmp.x();
				(*out_)[1] = tmp.y();
				(*out_)[2] = tmp.z();
			}
		}


		//-------------------------------------------------for Matrix
		///グローバル関数
		namespace function
		{

			///転置
			template<typename T, int Row, int Column>
			inline void MakeMatrixTranspose(
				const TMatrix<T, Row, Column>& matrix_,
				TMatrix<T, Column, Row>* out_)
			{
				assert(&matrix_ != out_);

				for(int c=0; c<Column; c++)
				for(int r=0; r<Row; r++)
					(*out_)(r, c) = matrix_(c, r);
			}
			
			///零行列
			template<typename T, int Row, int Column>
			inline void MakeMatrixZero(TMatrix<T, Row, Column>* out_)
			{
				const T zero(0);
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					(*out_)(r, c) = zero;
			}

			///単位行列(項等行列)
			template<typename T, int Row, int Column>
			inline void MakeMatrixIdentity(TMatrix<T, Row, Column>* out_)
			{
				const T unit(1);
				const T zero(0);
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					(*out_)(r, c)  = r == c ? unit : zero;
			}

			///拡大縮小
			template<template<typename, int> class Vector, typename T, int Row, int Column>
			inline void MakeMatrixScaling(
				const Vector<T, Column-1>& v_,
				TMatrix<T, Row, Column>* out_)
			{
				assert(Column == Row);
				TMatrix<T, Row, Column> obj;
				MakeMatrixIdentity(&obj);
				for(int i=1; i<Column-1; i++)
					(*out_)(i, i) = v_[i];
			}

			///回転X
			template<typename T>
			inline void MakeMatrixRotationX(const T& angle_, TMatrix<T, 4, 4>* out_)
			{
				float sin = function::Sin(angle_);
				float cos = function::Cos(angle_);
				
				MakeMatrixIdentity(out_);
				(*out_)(1, 1) =  cos;
				(*out_)(1, 2) = -sin;
				(*out_)(2, 1) =  sin;
				(*out_)(2, 2) =  cos;
			}		

			
			///回転Y
			template<typename T>
			inline void MakeMatrixRotationY(const T& angle_, TMatrix<T, 4, 4>* out_)
			{
				float sin = function::Sin(angle_);
				float cos = function::Cos(angle_);
				
				MakeMatrixIdentity(out_);
				(*out_)(0, 0) =  cos;
				(*out_)(0, 2) = -sin;
				(*out_)(2, 0) =  sin;
				(*out_)(2, 2) =  cos;
			}
				
			///回転Z
			template<typename T>
			inline void MakeMatrixRotationZ(const T& angle_, TMatrix<T, 4, 4>* out_)
			{
				float sin = function::Sin(angle_);
				float cos = function::Cos(angle_);
				
				MakeMatrixIdentity(out_);
				(*out_)(0, 0) =  cos;
				(*out_)(0, 1) = -sin;
				(*out_)(1, 0) =  sin;
				(*out_)(1, 1) =  cos;
			}

			///平行移動
			template<template<typename, int> class Vector, typename T, int MatrixDim>
			void MakeMatrixTranslation(
				const Vector<T, MatrixDim-1>& v_,
				TMatrix<T, MatrixDim, MatrixDim>* out_)
			{
				MakeMatrixIdentity(out_);
				for(int i=0; i<MatrixDim-1; i++)
					(*out_)(MatrixDim-1, i) = v_(i);
			}

			///せん断X
			template<typename T>
			inline void MakeMatrixShearX(const T& y_, const T& z_, const T& angle_, TMatrix<T, 4, 4>* out_)
			{
				MakeMatrixIdentity(out_);
				(*out_)(0, 1) = y_;
				(*out_)(0, 2) =	z_;

			}

			///せん断Y
			template<typename T>
			inline void MakeMatrixShearY(const T& z_, const T& x_, TMatrix<T, 4, 4>* out_)
			{
				MakeMatrixIdentity(out_);
				(*out_)(1, 1) = z_;
				(*out_)(1, 2) =	x_;

			}

			///せん断Z
			template<typename T>
			inline void MakeMatrixRotationZ(const T& x_, const T& y_, const T& angle_, TMatrix<T, 4, 4>* out_)
			{
				MakeMatrixIdentity(out_);
				(*out_)(2, 1) = y_;
				(*out_)(2, 2) =	z_;

			}

			///同時変換行列から平行移動成分を取り出す
			template<typename T>
			inline void GetTranslation(const TMatrix<T, 4, 4>& m_, TVector<T, 3>* out_)
			{
				(*out_)(0) = m_(3, 0);
				(*out_)(1) = m_(3, 1);
				(*out_)(2) = m_(3, 2);
			}

			///同時変換行列から回転成分のX方向ベクトルを取り出す
			template<typename T>
			inline void GetRotDirectionX(const TMatrix<T, 4, 4>& m_, TVector<T, 3>* out_)
			{
				(*out_)(0) = m_(0, 0);
				(*out_)(1) = m_(0, 1);
				(*out_)(2) = m_(0, 2);
			}

			///同時変換行列から回転成分のYX方向ベクトルを取り出す
			template<typename T>
			inline void GetRotDirectionY(const TMatrix<T, 4, 4>& m_, TVector<T, 3>* out_)
			{
				(*out_)(0) = m_(1, 0);
				(*out_)(1) = m_(1, 1);
				(*out_)(2) = m_(1, 2);
			}
			
			///同時変換行列から回転成分のZ方向ベクトルを取り出す
			template<typename T>
			inline void GetRotDirectionZ(const TMatrix<T, 4, 4>& m_, TVector<T, 3>* out_)
			{
				(*out_)(0) = m_(2, 0);
				(*out_)(1) = m_(2, 1);
				(*out_)(2) = m_(2, 2);
			}

			///3×3部分の行列式を得る
			template<typename T>
			inline void Determinant33(const TMatrix<T, 4, 4>& m_, T* out_)
			{
				*out_ = m_(0, 0)*(m_(1, 1)*m_(2, 2) - m_(1, 2)*m_(2, 1))
				      + m_(0, 1)*(m_(1, 2)*m_(2, 0) - m_(1, 0)*m_(2, 2))
				      + m_(0, 2)*(m_(1, 0)*m_(2, 1) - m_(1, 1)*m_(2, 0));
			}

			///4×3部分を逆行列にする                     <br>
			///<br>
			///逆行列の公式により直接求める（掃出法より高速と判断）<br>
			///<br>
			///Inv(M) = adj(M)/det(M)                          <br>
			///<br>
			template<typename T>
			inline void Inverse43(const TMatrix<T, 4, 4>& m_, TMatrix<T, 4, 4>* out_)
			{
				//行列式
				T det;
				Determinant33(m_, det);

				//特異行列の0割チェック
				assert(function::Abs(det) > 0.000001f);

				//行列式の逆数を得る
				T overDet = 1 / det;

				//逆行列はadj(M)/det(M)
				(*out_)(0, 0) = (m_(1, 1)*m_(2, 2) - m_(1, 2)*m_(2, 1))*overDet;
				(*out_)(0, 1) = (m_(0, 2)*m_(2, 1) - m_(0, 1)*m_(2, 2))*overDet;
				(*out_)(0, 2) = (m_(0, 1)*m_(1, 2) - m_(0, 2)*m_(1, 1))*overDet;

				(*out_)(1, 0) = (m_(1, 2)*m_(2, 0) - m_(1, 0)*m_(2, 2))*overDet;
				(*out_)(1, 1) = (m_(0, 0)*m_(2, 2) - m_(0, 2)*m_(2, 0))*overDet;
				(*out_)(1, 2) = (m_(0, 2)*m_(1, 0) - m_(0, 0)*m_(1, 2))*overDet;

				(*out_)(2, 0) = (m_(1, 0)*m_(2, 1) - m_(1, 1)*m_(2, 0))*overDet;
				(*out_)(2, 1) = (m_(0, 1)*m_(2, 0) - m_(0, 0)*m_(2, 1))*overDet;
				(*out_)(2, 2) = (m_(0, 0)*m_(1, 1) - m_(0, 1)*m_(1, 0))*overDet;

				(*out_)(3, 0) = -(m_(3, 0)*m_(0, 0) + m_(3, 1)*m_(1, 0) + m_(3, 2)*m_(2, 0));
				(*out_)(3, 1) = -(m_(3, 0)*m_(0, 1) + m_(3, 1)*m_(1, 1) + m_(3, 2)*m_(2, 1));
				(*out_)(3, 2) = -(m_(3, 0)*m_(0, 2) + m_(3, 1)*m_(1, 2) + m_(3, 2)*m_(2, 2));
			}

			///回転軸と回転角度から回転行列を生成                                   <br>
			///<br>
			///n ：回転軸方向の単位ベクトル                                                 <br>
			///v ：回転前                                                                   <br>
			///v'：回転後                                                                   <br>
			///θ：回転角                                                                   <br>
			///<br>
			///vを回転軸と垂直な方向v'⊥と平行な方向v'∥に分解し                            <br>
			///<br>
			///v = v⊥ + v∥                                                            <br>
			///<br>
			///とすると                                                                     <br>
			///<br>
			///v∥ = (n・v)*n                                                           <br>
			///v⊥ = v - v∥                                                            <br>
			///<br>
			///vとnに垂直なベクトルwとすると                                                <br>
			///<br>
			///w = n × v                                                               <br>
			///<br>
			///このとき回転後のベクトルについて垂直な方向v'⊥と平行な方向v'∥とすると       <br>
			///<br>
			///v'⊥ = v⊥*cosθ + w*sinθ                                               <br>
			///v'∥ = v∥                                                               <br>
			///<br>
			///よって回転後のベクトルは v' = v'⊥ + v'∥ により                             <br>
			///<br>
			///v' = (v - (n・v)*n)*cosθ + (n × v)*sinθ + (n・v)n                     <br>
			///<br>
			///この変換をex=[1 0 0], ey=[0 1 0], ez=[0 0 1]のそれぞれに対して行い           <br>
			///行ベクトルを並べることで任意軸周りの回転行列を得る                           <br>
			///<br>
			///ex' = [nx^2+(1-nx^2)*cos,    nx*ny*(1-cos)+z*sin,  nx*ny*(1-cos)-ny*sin] <br>
			///ey' = [nx*ny*(1-cos)-nz*sin, ny^2+(1-ny^2)*cos,    ny*nz*(1-cos)+nx*sin] <br>
			///ez' = [nx*nz*(1-cos)+ny*sin, ny*nz*(1-cos)-nx*sin, nz^2+(1-nz^2)*cos   ] <br>
			///<br>
			///∴M = [ex';                                                              <br>
			///ey';                                                              <br>
			///ez']                                                              <br>
			///<br>
			template<template<typename, int>class Vector, typename T>
			inline void MakeMatrixFromAxisAngle(
				const Vector<T, 3>& axis_,
				const T&  angle_,
				TMatrix<T, 4, 4>* out_)
			{
				if (function::Abs(axisIdentity_.LengthSquared() - 1.0) > 0.01)
				{
					assert("axis must be normalized\n");
				}

				float sin = Sin(angle_);
				float cos = Cos(angle_);
	             
				float i = 1 - cos;
				float xs = axis_[0] * sin;
				float ys = axis_[1] * sin;
				float zs = axis_[2] * sin;
				float xi = axis_[0] * i;
				float yi = axis_[1] * i;
				float zi = axis_[2] * i;
				float xx = xi * axis_[0];
				float yy = yi * axis_[1];
				float zz = zi * axis_[2];
				float xy = xi * axis_[1];
				float yz = yi * axis_[2];
				float zx = zi * axis_[0];
			
				(*out_)(0, 0) = xx + cos;
				(*out_)(0, 1) = xy + zs;
				(*out_)(0, 2) = zx - ys;
				(*out_)(0, 3) = 0;
				(*out_)(1, 0) = xy - zs;
				(*out_)(1, 1) = yy + cos;
				(*out_)(1, 2) = yz + xs;
				(*out_)(1, 3) = 0;
				(*out_)(2, 0) = zx + ys;
				(*out_)(2, 1) = yz - xs;
				(*out_)(2, 2) = zz + cos;
				(*out_)(2, 3) = 0;
				(*out_)(3, 0) = 0;
				(*out_)(3, 1) = 0;
				(*out_)(3, 2) = 0;
				(*out_)(3, 3) = 1;
			}

			///任意軸方向へのスケーリング行列を生成 　                    <br>
			///<br>
			///n ：軸方向の単位ベクトル                                           <br>
			///v ：スケーリング前                                                 <br>
			///v'：スケーリング後                                                 <br>
			///k ：係数                                                           <br>
			///<br>
			///vを軸と垂直な方向v'⊥と平行な方向v'∥に分解し                      <br>
			///<br>
			///v = v⊥ + v∥                                                  <br>
			///<br>
			///とすると                                                           <br>
			///<br>
			///v∥ = (n・v)n                                                  <br>
			///v⊥ = v - v∥                                                  <br>
			///<br>
			///このときスケーリング後のベクトルについて                           <br>
			///垂直な方向v'⊥と平行な方向v'∥とすると                             <br>
			///<br>
			///v'⊥ = v⊥                                                     <br>
			///v'∥ = k*v∥                                                   <br>
			///<br>
			///よってスケーリング後のベクトルは v' = v'⊥ + v'∥ により           <br>
			///<br>
			///v' = v - (n・v)*n + k*(n・v)n                                  <br>
			///<br>
			///この変換をex=[1 0 0], ey=[0 1 0], ez=[0 0 1]のそれぞれに対して行い <br>
			///行ベクトルを並べることで任意軸方向のスケーリング行列を得る         <br>
			///<br>
			///ex' = [1+(k-1)*nx^2, (k-1)*nx*ny,  (k-1)*nx*nz ]               <br>
			///ey' = [(k-1)*nx*ny,  1+(k-1)*ny^2, (k-1)*ny*nz ]               <br>
			///ez' = [(k-1)*nx*nz,  (k-1)*ny*nx,  1+(k-1)*nz^2]               <br>
			///<br>
			///∴M = [ex';                                                    <br>
			///ey';                                                    <br>
			///ez']                                                    <br>
			///<br>
			template<template<typename, int>class Vector, typename T>              
			inline void MakeMatrixFromAxisScale(                                   
				const Vector<T, 3>& axisIdentity_,
				const T&  k_,
				TMatrix<T, 4, 4>* out_)
			{
				if (function::Abs(axisIdentity_.LengthSquared() - 1.0) > 0.01)
				{
					assert("axis must be normalized\n");
				}

				float k1  = k_ - 1;
				float kxx = 1 + k1 * axisIdentity_[0] * axisIdentity_[1];
				float kyy = 1 + k1 * axisIdentity_[1] * axisIdentity_[2];
				float kzz = 1 + k1 * axisIdentity_[2] * axisIdentity_[0];
				float kxy = k1 * axisIdentity_[0] * axisIdentity_[1];
				float kyz = k1 * axisIdentity_[1] * axisIdentity_[2];
				float kzx = k1 * axisIdentity_[2] * axisIdentity_[0];

				(*out_)(0, 0) = kxx;
				(*out_)(0, 1) = kxy;
				(*out_)(0, 2) = kzx;
				(*out_)(0, 3) = 0;
				(*out_)(1, 0) = kxy;
				(*out_)(1, 1) = kyy;
				(*out_)(1, 2) = kyz;
				(*out_)(1, 3) = 0;
				(*out_)(2, 0) = kzx;
				(*out_)(2, 1) = kyz;
				(*out_)(2, 2) = kzz;
				(*out_)(2, 3) = 0;
				(*out_)(3, 0) = 0;
				(*out_)(3, 1) = 0;
				(*out_)(3, 2) = 0;
				(*out_)(3, 3) = 1;
			}

						
			///オイラー角から回転を含んだ4×4同次変換行列に変換                  <br>
			///<br>
			///平行移動の成分は                                                          <br>
			///pitch→heading→ybankをx軸→x軸→z軸の順番に適用した行列の結合と同じ      <br>
			///<br>
			///[cosh*cosb + sinh*sinp*sinb,  sinb*cosp, -sinh*cosb + cosh*sinp*sinb; <br>
			///-cosh*sinb + sinh*sinp*cosb, cosb*cosp, sinb*sinh + cosh*sinp*cosb ; <br>
			///sinh*cisp,                   -sinp,     cosh*cosp                  ] <br>
			///<br>
			template<typename T>
			inline void MakeMatrixFromEulerAngles(
				const TEulerAngles<T>& e_,
				TMatrix<T, 4, 4>* out_)
			{
				T cosp = function::Cos(e_.pitch());
				T cosh = function::Cos(e_.heading());
				T cosb = function::Cos(e_.bank());
				T sinp = function::Sin(e_.pitch());
				T sinh = function::Sin(e_.heading());
				T sinb = function::Sin(e_.bank());
				
				(*out_)(0, 0) = cosh*cosb + sinh*sinp*sinb;
				(*out_)(0, 1) = sinb*cosp;
				(*out_)(0, 2) = -sinh*cosb + cosh*sinp*sinb;
				(*out_)(0, 3) = T(0);
				(*out_)(1, 0) = -cosh*sinb + sinh*sinp*cosb;
				(*out_)(1, 1) = cosb*cosp;
				(*out_)(1, 2) = sinb*sinh + cosh*sinp*cosb;
				(*out_)(1, 3) = T(0);
				(*out_)(2, 0) = sinh*cosp;
				(*out_)(2, 1) = -sinp;
				(*out_)(2, 2) = cosh*cosp;
				(*out_)(2, 3) = T(0);
				(*out_)(3, 0) = T(0);
				(*out_)(3, 1) = T(0);
				(*out_)(3, 2) = T(0);
				(*out_)(3, 3) = T(1);
			}

			///クォータニオンから回転を含んだ4×4同次変換行列に変換               <br>
			///<br>
			///任意軸周りの回転行列                                                       <br>
			///<br>
			///R = [nx^2+(1-nx^2)*cos,    nx*ny*(1-cos)+nz*sin, nx*ny*(1-cos)-ny*sin; <br>
			///nx*ny*(1-cos)-nz*sin, ny^2+(1-ny^2)*cos,    ny*nz*(1-cos)+nx*sin; <br>
			///nx*nz*(1-cos)+ny*sin, ny*nz*(1-cos)-nx*sin, nz^2+(1-nz^2)*cos   ] <br>
			///<br>
			///任意軸周りの回転クォータニオン                                             <br>
			///<br>
			///x = nx*sin(θ/2)                                                       <br>
			///y = ny*sin(θ/2)                                                       <br>
			///z = nz*sin(θ/2)                                                       <br>
			///w = cos(θ/2)                                                          <br>
			///<br>
			///これらを比較して回転クォータニオンを逆算する。                             <br>
			///対角要素はcos倍角公式 cos(2*θ) = 1- 2sin(θ)^2 により                     <br>
			///<br>
			///R11 = 1 - (1-nx^2)(1-cosθ)                                            <br>
			///= 1 - (1-nx^2)(2*sin(θ/2)^2)                                      <br>
			///= -1 + 2*cos(θ/2)^2 + 2*(nx*cos(θ/2))^2                          <br>
			///= -1 + 2*w^2 + 2*x^2                                               <br>
			///<br>
			///またx^2 + y^2 + z^2 + w^2 = 1により以下も同値である                        <br>
			///<br>
			///R11 = 1-(1-nx^2)(1-cosθ)                                              <br>
			///= 1-(1-nx^2)(2*sin(θ/2)^2)                                        <br>
			///= 1-(ny^2+nz^2)(2*sin(θ/2)^2)   ※                                 <br>
			///= 1-2*cos(θ/2)^2 - 2*(nx*cos(θ/2))^2                             <br>
			///= 1-2*y^2 - 2*z^2                                                  <br>
			///<br>
			///非対角要素はsin倍角公式 sin(2*θ) = 1- 2sinθcosθにより                   <br>
			///<br>
			///R12 = nx*ny*(1-cosθ)+nz*sinθ                                         <br>
			///= nx*ny*(2*sin(θ/2)^2)+nz*2*sin(θ/2)cos(θ/2)                    <br>
			///= 2xy + 2wz                                                        <br>
			///<br>
			///以下同様にして                                                             <br>
			///R = [1 - 2*y^2 - 2*z^2, 2xy + 2wz,         2xz - 2wy        ;          <br>
			///2xy - 2wz,         1 - 2*x^2 - 2*z^2, 2yz + 2wx        ;          <br>
			///2xz + 2wy,         2yz - 2wx,         1 - 2*x^2 - 2*y^2]          <br>
			///<br>
			template<typename T>
			inline void MakeMatrixFromQuaternion(
				const TQuaternion<T>& q_, 
				TMatrix<T, 4, 4>* out_)
			{
				TQuaternion<T> qw(
					q_.x()*q_.x(),
					q_.y()*q_.y(),
					q_.z()*q_.z(),
					q_.w()*q_.w());
				qw += qw;

				TQuaternion<T> q2(q_);
				q2 += q2;

				T xy2 = q2.x() * q_.y();
				T yz2 = q2.y() * q_.z();
				T zx2 = q2.x() * q_.z();
				T xw2 = q_.x() * q2.w();
				T yw2 = q_.y() * q2.w();
				T zw2 = q_.z() * q2.w();
				
				(*out_)(0, 0) = 1 - qw.y() - qw.z();
				(*out_)(0, 1) = xy2 + zw2;
				(*out_)(0, 2) = zx2 - yw2;
				(*out_)(0, 3) = 0;
				(*out_)(1, 0) = xy2 - zw2;
				(*out_)(1, 1) = 1 - qw.z() - qw.x();
				(*out_)(1, 2) = yz2 + xw2;
				(*out_)(1, 3) = 0;
				(*out_)(2, 0) = zx2 + yw2;
				(*out_)(2, 1) = yz2 - xw2;
				(*out_)(2, 2) = 1 - qw.x() - qw.y();
				(*out_)(2, 3) = 0;
				(*out_)(3, 0) = 0;
				(*out_)(3, 1) = 0;
				(*out_)(3, 2) = 0;
				(*out_)(3, 3) = 1;
			}
		}


		//-------------------------------------------------for Euler Angle
		///グローバル関数
		namespace function
		{
			///項等オイラー角
			template<typename T>
			inline void MakeEulerAnglesIdentity(TEulerAngles<T>* out_)
			{
				const T zero(0);
				out_->pitch()   = zero(0);
				out_->heading() = zero(0);
				out_->bank()    = zero(0);
			}

			///正準化
			template<typename T>
			inline void MakeEulerAnglesCanonize(
				const TEulerAngles<T>& m_,
				TEulerAngles<T>* out_)
			{
				(*out_) = m_;

				//pitchを正準範囲[-π/2, π/2]に変換する
				//[-π, π]の範囲に周回した後に裏表を判定して正準範囲に変換
				math::function::Cycle(&out_->pitch(), -definition::RF_PI, definition::RF_PI); 
				if (out_->pitch() < -definition::RF_OVER_2PI)
				{
					out_->pitch()   -= definition::RF_PI;
					out_->heading() += definition::RF_PI;
					out_->bank()    += definition::RF_PI;
				}
				else if(definition::OVER_PI2 < out_->pitch)
				{
					out_->pitch()    = definition::RF_PI - out_->pitch();
					out_->heading() += definition::RF_PI;
					out_->bank()    += definition::RF_PI;
				}

				//gimbal lockしていなければ
				//headingを正準範囲に周回
				if (std::fabs(out_->pitch())>definition::RF_OVER_2PI - definition::RF_EPSILON)
				{
					out_->heading() += m_.bank();
					out_->bank() = 0.0f;
				}
				else
				{
					math::function::Cycle(&out_->bank(), -definition::RF_PI, definition::RF_PI); 
				}

				//headingを正準範囲に周回
				math::function::Cycle(&out_->heading(), -definition::RF_PI, definition::RF_PI); 

				return *this;
			}

			
			///回転を含んだ4×4同次変換行列からオイラー角に変換                      <br>
			///<br>
			///オイラー角から行列への逆変換                                                  <br>
			///<br>
			///R = [cosh*cosb + sinh*sinp*sinb,  sinb*cosp, -sinh*cosb + cosh*sinp*sinb; <br>
			///-cosh*sinb + sinh*sinp*cosb, cosb*cosp, sinb*sinh + cosh*sinp*cosb ; <br>
			///sinh*cisp,                   -sinp,     cosh*cosp                  ] <br>
			///<br>
			///まずpがすぐに求まる                                                           <br>
			///<br>
			///p = asin(-R32)                                                            <br>
			///<br>
			///次に                                                                        <br>
			///<br>
			///R31/cosp = sinh                                                           <br>
			///R33/cosp = cosh                                                           <br>
			///<br>
			///になるのでtan = sin/cosを利用して                                             <br>
			///<br>
			///h = atan(R31/cosp, R33/cosp) = atan(R31, R33)                             <br>
			///<br>
			///同様に                                                                        <br>
			///<br>
			///b = atan(R12, R22)                                                        <br>
			///<br>
			///ただしジンバルロックしている場合を考慮する場合                                <br>
			///<br>
			///cosp ≒ 0なら b = 0とする                                                 <br>
			///<br>
			///としてそれぞれRに代入すると                                                   <br>
			///<br>
			///R = [cosh,      0,     -sinh     ;                                        <br>
			///sinh*sinp, 0,     cosh*sinp ;                                        <br>
			///0,         -sinp, 0         ]                                        <br>
			///<br>
			///となり、R11とR13からこのケースの回転角度を求めることとする                    <br>
			///<br>
			template<typename T>
			inline void MakeEulerAnglesFromMatrix(
				const TMatrix<T, 4, 4>& m_,
				TEulerAngles<T>* out_)
			{
				T sp = -m_(2, 1);

				//sinpが±1付近→pitchが±PI/2→ジンバルロック
				if (function::Abs(sp) > T(1) - T(definition::RF_EPSILON))
				{
					out_->pitch()   = definition::RF_PI_OVER2 * sp; //?
					out_->heading() = function::ATan2(-m_(1, 2), m_(0, 0));
					out_->bank()    = 0.0f;
				}
				else
				{
					out_->pitch()   = function::ASin(sp);
					out_->heading() = function::ATan2(m_(2, 0), m_(2, 2));
					out_->bank()    = function::ATan2(m_(0, 1), m_(1, 1));
				}
			}

			///クォータニオンをオイラー角に変換                          <br>
			///<br>
			///行列からオイラー角への変換                                        <br>
			///<br>
			///if (gimbalLock){                                              <br>
			///pitch   = PI/2 * (-m32);                                  <br>
			///heading = atan2(-m23, m11);                               <br>
			///bank    = 0.0f;                                           <br>
			///}else{                                                        <br>
			///pitch   = ASin(-m32);                                     <br>
			///heading = ATan2(m31, m_33);                               <br>
			///bank    = ATan2(m12, m22);                                <br>
			///}                                                             <br>
			///<br>
			///に以下を代入して求めることが可能。                                <br>
			///<br>
			///R = [1 - 2*y^2 - 2*z^2, 2xy + 2wz,         2xz - 2wy        ; <br>
			///2xy - 2wz,         1 - 2*x^2 - 2*z^2, 2yz + 2wx        ; <br>
			///2xz + 2wy,         2yz - 2wx,         1 - 2*x^2 - 2*y^2] <br>
			///<br>
			///ちなみに、以下オイラー角からクォータニオンを求める                <br>
			///式から解くことも可能らしい。                                      <br>
			///<br>
			///w =  cosh*cosp*cosb + cosh*sinp*sinb;                         <br>
			///x =  cosh*sinp*cosb + sinh*cosp*sinb;                         <br>
			///y = -cosh*sinp*sinb + sinh*cosp*cosb;                         <br>
			///z = -cosh*sinp*cosb + cosh*cosp*sinb;                         <br>
			///<br>
			template<typename T>
			inline void MakeEulerAnglesFromQuaternion(
				const TQuaternion<T>& q_,
				TEulerAngles<T>* out_)
			{
				T sp = -2.0f * (q_.y()*q_.z() - q_.w()*q_.x());

				if(function::Abs(sp) > T(1) - T(definition::RF_EPSILON))
				{
					out_->pitch() = definition::RF_OVER_2PI * sp;

					out_->heading() = function::ATan2(
						-q_.x()*q_.z() + q_.w()*q_.y(),
						0.5f - q_.y()*q_.y() - q_.z()*q_.z());

					out_->bank() = 0.0f;
				}
				else
				{
					out_->pitch() = function::ASin(sp);

					out_->heading() = function::ATan2(
						-q_.x()*q_.z() + q_.w()*q_.y(),
						0.5f - q_.x()*q_.x() - q_.y()*q_.y());

					out_->bank() = function::ATan2(
						-q_.x()*q_.y() + q_.w()*q_.z(),
						0.5f - q_.x()*q_.x() - q_.z()*q_.z());
				}
			}
		}


		//-------------------------------------------------for Quaternion
		///グローバル関数
		namespace function
		{
			///内積
			template<typename T>
			T Dot(const TQuaternion<T>& lhs_, const TQuaternion<T>& rhs_)
			{
				return lhs_.x()*rhs_.x() + lhs_.y()*rhs_.y() + lhs_.z()*rhs_.z() + lhs_.w()*rhs_.w();
			}

			///累乗                               <br>
			///<br>
			///クォータニオンは累乗すると回転が加算される <br>
			///<br>
			///exp(log(q)) = q                            <br>
			///q^t = exp(t*log(q)) ..etc..                <br>
			///<br>
			///q^(st) ≠ q^s^t                            <br>
			///<br>
			template<typename T>
			inline void MakeQuaternionPow(const TQuaternion<float>&src_, int exp_, TQuaternion<float>* out_)
			{
				if (function::Abs(src_.w()) > 0.9999f)
				{
					LOCALPRINTF("quaternion identity? %s\n", ((std::string)src_).c_str());
					return;
				}

				//θを取得し、指数を掛ける
				T theta = function::ACos(src_.w());
				T newTheta = theta * exp_;

				T mult = function::Sin(newTheta) / function::Sin(theta)

				out_->x() = src_.x() * mult;
				out_->y() = src_.y() * mult;
				out_->z() = src_.z() * mult;
				out_->w() = function::Cos(newTheta);
			}

			///pure quaternion
			template<typename T>
			inline void MakeQuaternionPure(const TVector<T, 3>&src_, TQuaternion<T>* dest_)
			{
				dest_->x() = src_(0);
				dest_->y() = src_(1);
				dest_->z() = src_(2);
				dest_->w() = T(0);
			}

			///項等クォータニオン
			template<typename T>
			inline void MakeQuaternionIdentity(TQuaternion<T>* dest_)
			{
				dest_->x() = T(0);
				dest_->y() = T(0);
				dest_->z() = T(0);
				dest_->w() = T(1);
			}

			///共役クォータニオン
			template<typename T>
			inline void MakeQuaternionConjugate(const TQuaternion<T>& src_, TQuaternion<T>* dest_)
			{
				dest_->x() = -src_.x();
				dest_->y() = -src_.y();
				dest_->z() = -src_.z();
				dest_->w() = src_.w();
			}
			
			///逆クォータニオン
			template<typename T>
			inline void MakeQuaternionInverse(const TQuaternion<T>& src_, TQuaternion<T>* dest_)
			{
				T ls = src_.ReciprocalLength();
				ls *= ls;
				dest_->x() = -src_.x()*ls;
				dest_->y() = -src_.y()*ls;
				dest_->z() = -src_.z()*ls;
				dest_->w() = src_.w()*ls;
			}
			
			///正規化
			template<typename T>
			inline void MakeQuaternionNormalize(const TQuaternion<T>& src_, TQuaternion<T>* dest_)
			{
				T rl = src_.ReciprocalLength();
				dest_->x() = src_.x()*rl;
				dest_->y() = src_.y()*rl;
				dest_->z() = src_.z()*rl;
				dest_->w() = src_.w()*rl;
			}
			
			///回転X
			template<typename T>
			inline void MakeQuaternionRotationX(const T& angle_, TQuaternion<T>* out_)
			{
				T tmp = angle_ * 0.5f;

				out_->x() = function::Sin(tmp);
				out_->y() = 0;
				out_->z() = 0;
				out_->w() = function::Cos(tmp);
			}

			///回転Y
			template<typename T>
			inline void MakeQuaternionRotationY(const T& angle_, TQuaternion<T>* out_)
			{
				T tmp = angle_ * 0.5f;

				out_->x() = 0;
				out_->y() = function::Sin(tmp);
				out_->z() = 0;
				out_->w() = function::Cos(tmp);
			}
				
			///回転Z
			template<typename T>
			inline void MakeQuaternionRotationZ(const T& angle_, TQuaternion<T>* out_)
			{
				T tmp = angle_ * 0.5f;

				out_->x() = 0;
				out_->y() = 0;
				out_->z() = function::Sin(tmp);
				out_->w() = function::Cos(tmp);
			}

			///任意軸周りの回転                        <br>
			///使う側はaxis_を正規化しておく必要がある <br>
			template<typename T>
			inline void MakeQuaternionFromAxisAngle(
				const TVector<T, 3>& axisIdentity_, T angle_, TQuaternion<T>* out_)
			{
				//正規化チェック
				assert(function::Abs(axisIdentity_.LengthSquared() - T(1.0)) < T(0.01));

				angle_ *= 0.5f;
				T sin = Sin(angle_);
				T cos = Cos(angle_);
				out_->x() = -axisIdentity_[0]*sin;
				out_->y() = -axisIdentity_[1]*sin;
				out_->z() = -axisIdentity_[2]*sin;
				out_->w() = cos;
			}

			///Linear Interpolation
			template<typename T>
			inline void MakeQuaternionLerp(
				const TQuaternion<T>& a_,
				const TQuaternion<T>& b_,
				const typename TQuaternion<T>::value_type& amount_,
				TQuaternion<T>* out_)
			{
				*out_ = a_*(1 - amount_) + b_*amount_;
			}

			///Slerp (Spherical Linear IntERPolation)                                <br>  
			///<br>
			///クォータニオンは累乗することで回転を重ねることが出来るため                    <br>
			///q0からq1への角変移は以下                                                      <br>
			///<br>
			///Δq = Inv(q0) * q1                                                        <br>
			///<br>
			///このフラクション差分はΔq^t                                                   <br>
			///よってSlerp(q0, q1, t) = q0 * (Inv(q0) * q1)^t                                <br>
			///が、この計算は処理が重いので、                                                <br>
			///ベクトルの補間法を使う                                                        <br>
			///<br>
			///vt     = k0*v0 + k1*v1                                                    <br>
			///<br>
			///v1                                                                   <br>
			////    k1*v1                                                            <br>
			////    /                                                                 <br>
			////    /                                                                  <br>
			////    /                                                                   <br>
			///o/____/______ v0                                                           <br>
			///k0v0                                                                     <br>
			///<br>
			///とし、v0とv1がなす角θとし、                                                  <br>
			///v0とk1*v1がなす角はt*θとすると、k1は以下                                     <br>
			///<br>
			///sin(θ) / 1 = sin(t*θ) / k1により                                        <br>
			///k1 = sin(t*θ) / sin(θ)                                                  <br>
			///<br>
			///同様にしてk0は以下                                                            <br>
			///<br>
			///sin(θ) / 1 = sin(1-t) / k0                                               <br>
			///k0 = sin((1-t)*θ) / sin(θ)                                              <br>
			///<br>
			///∴                                                                            <br>
			///ベクトルの場合                                                            <br>
			///vt(v0, v1, t) = (sin((1-t)*θ)/sin(θ))*v0 + ((sin(t*θ) / sin(θ))*v1    <br>
			///<br>
			///４元数の場合も同様                                                        <br>
			///slerp(q0, q1, t) = (sin((1-t)*θ)/sin(θ))*q0 + ((sin(t*θ) / sin(θ))*q1 <br>
			///<br>
			///ただしθ≒0では　分母が発散するので線形補間とする。                           <br>
			///sin((1-t)θ)*a + sin(tθ)*b)                                              <br>
			///<br>
			template<typename T>
			inline void MakeQuaternionSlerp(
				const TQuaternion<T>& a_,
				const TQuaternion<T>& b_,
				const typename TQuaternion<T>::value_type& amount_,
				TQuaternion<T>* out_)
			{
				//正規化チェック
				assert(function::Abs(a_.LengthSquared() - 1.0) < 0.01);
				assert(function::Abs(b_.LengthSquared() - 1.0) < 0.01);

				const T RF_EPSILON(1e-4f);

				T cos = Dot(a_, b_);
				
				//原点付近では分母0になり計算が発散するので線形補間する
				if (cos >= 1 - RF_EPSILON)
				{
					MakeQuaternionLerp(a_, b_, amount_, out_);
				}
				else
				{
					if (-1 > cos)
						cos = -1.0f;
					else if(1 < cos)
						cos = 1.0f;

					T theta = amount_ * function::ACos(cos);
					T u = Cos(theta);
					T v = Sin(theta);
					TQuaternion<T> q2(b_ - a_*cos);
					q2.Normalize();
					q2 *= v;
					TQuaternion<T> q0(a_*u);
					*out_ = q0 + q2;
				}
			}

			///回転を含んだ4×4同次変換行列からクォータニオンに変換      <br>
			///<br>
			///R = [1 - 2*y^2 - 2*z^2, 2xy + 2wz,         2xz - 2wy        ; <br>
			///2xy - 2wz,         1 - 2*x^2 - 2*z^2, 2yz + 2wx        ; <br>
			///2xz + 2wy,         2yz - 2wx,         1 - 2*x^2 - 2*y^2] <br>
			///<br>
			///１）この対角要素の和（トレース）をとり                            <br>
			///w^2 + x^2 + 2*y^2 + 2*z^2 = 1を適用すると                     <br>
			///<br>
			///R11 + R22 + R33  = 3 - 4*(1 - w)                              <br>
			///R11 - R22 - R33  = 4*(x - 1)                                  <br>
			///-R11 + R22 - R33 = 4*(y - 1)                                  <br>
			///-R11 - R22 - R33 = 4*(z - 1)                                  <br>
			///<br>
			///ただし、これだけでは符号が判別できない。                          <br>
			///<br>
			///そこで、四元数は幾何学的には q=-q という性質を利用する。          <br>
			///即ち、任意の要素ひとつだけについては符号を無視することが出来る。  <br>
			///<br>
			///２）対角について対称の位置の和と差をとると                        <br>
			///<br>
			///R12 + R21 = 4xy                                               <br>
			///R12 - R21 = 4wz                                               <br>
			///R23 + R32 = 4yz                                               <br>
			///R23 - R32 = 4wx                                               <br>
			///R31 + R13 = 4zy                                               <br>
			///R31 - R13 = 4wy                                               <br>
			///<br>
			///１）の最大値によって基準となる要素を選択し値をひとつ決める。      <br>
			///２）によって符号を考慮した要素の値を得る                          <br>
			///<br>
			template<typename T>
			inline void MakeQuaternionFromMatrix(
				const TMatrix<T, 4, 4>& m_,
				TQuaternion<T>* out_)
			{

				T ww = m_(0, 0) + m_(1, 1) + m_(2, 2);
				T wx = m_(0, 0) + m_(1, 1) + m_(2, 2);
				T wy = m_(1, 1) + m_(0, 0) + m_(2, 2);
				T wz = m_(2, 2) + m_(0, 0) + m_(1, 1);
				
				T biggest = ww;
				int ix = 0;
				if (wx > biggest)
				{
					biggest = wx;
					ix = 1;
				}
				else if (wy > biggest)
				{
					biggest = wy;
					ix = 2;				
				}
				else if (wz > biggest)
				{
					biggest = wz;
					ix = 3;
				}
				
				T biggestVal = function::Sqrt(biggest + 1) * 0.5f;
				T mult = 1.0f / (4.0f*biggestVal);
				
				switch(ix)
				{
					case 0:
						out_->x() = (m_(1, 2) - m_(2, 1)) * mult;
						out_->y() = (m_(2, 0) - m_(0, 2)) * mult;
						out_->z() = (m_(0, 1) - m_(1, 0)) * mult;
						out_->w() = biggestVal;
						break;
					case 1:
						out_->x() = biggestVal;
						out_->y() = (m_(0, 1) + m_(1, 0)) * mult;
						out_->z() = (m_(2, 0) + m_(0, 2)) * mult;
						out_->w() = (m_(1, 2) - m_(2, 1)) * mult;
						break;
					case 2:
						out_->x() = (m_(0, 1) + m_(1, 0)) * mult;
						out_->y() = biggestVal;
						out_->z() = (m_(1, 2) + m_(2, 1)) * mult;
						out_->w() = (m_(2, 0) - m_(0, 2)) * mult;
						break;
					case 3:
						out_->x() = (m_(2, 0) + m_(0, 2)) * mult;
						out_->y() = (m_(1, 2) + m_(2, 1)) * mult;
						out_->z() = biggestVal;
						out_->w() = (m_(0, 1) - m_(1, 0)) * mult;
						break;
				}
			}


			///オイラー角からクォータニオン          <br>
			///hpbの順にクォータニオンを掛ければよい <br>
			template<typename T>
			inline void MakeQuaternionFromEulerAngles(
				const TEulerAngles<T>& e_,
				TQuaternion<T>* out_)
			{
				T sinp;
				T sinh;
				T sinb;
				T cosp;
				T cosh;
				T cosb;

				SinCos(e_.pitch()*T(0.5),   &sinp, &cosp);
				SinCos(e_.heading()*T(0.5), &sinh, &cosh);
				SinCos(e_.bank()*T(0.5),    &sinb, &cosb);
		
				out_->w() =  cosh*cosp*cosb + cosh*sinp*sinb;
				out_->x() =  cosh*sinp*cosb + sinh*cosp*sinb;
				out_->y() = -cosh*sinp*sinb + sinh*cosp*cosb;
				out_->z() = -cosh*sinp*cosb + cosh*cosp*sinb;
			}
		}
	}
}

#endif

