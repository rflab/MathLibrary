//2010/04/18
//rflab.

#ifndef _RF_QUATERNION_TEMPLATE_H_
#define _RF_QUATERNION_TEMPLATE_H_

#include <cassert>
#include <sstream>

#include "math.h"
#include "expressiontemplate.h"

//debug
#include <stdio.h>
#undef LOCALPRINTF
//#define LOCALPRINTF(...) printf(__VA_ARGS__)
#define LOCALPRINTF(...)

///rock field!
namespace rf
{
	///数学
	namespace math
	{
		
		///クォータニオン・四元数<br>
		///<br>
		///ijk = i^2 = j^2 = k^2 = -1    <br>
		///ij = -ji = k                  <br>
		///jk = -kj = i                  <br>
		///<br>
		///m_w + (m_x*i + m_y*j + m_z*k) <br>
		///<br>
		///これを使っていろいろ出来る！  <br>
		///<br>
		template<typename T>
		class TQuaternion
		{
		public:
			T m_x; //imaginary i
			T m_y; //imaginary j
			T m_z; //imaginary k
			T m_w; //real
			
		public:
			typedef T value_type;
			typedef TQuaternion<T> self_type;
		public:
			
			///デストラクタ
			///virtualはつけない->速度が異様に落ちる
			~TQuaternion(){}

			///引数無しコンストラクタ
			///初期化も含めて何もしない
			TQuaternion()
			{
				LOCALPRINTF("TQuaternion no arg  constructor\n");
			}

			///値指定コンストラクタ
			TQuaternion(T x_, T y_, T z_, T w_)
			{
				LOCALPRINTF("TQuaternion each member init constructor\n");
				m_x = x_;
				m_y = y_;
				m_z = z_;
				m_w = w_;
			}

			///コピーコンストラクタ
			TQuaternion(const self_type& rhs_) 
			{
				LOCALPRINTF("TQuaternion copy constructor\n");
				m_x = rhs_.m_x;
				m_y = rhs_.m_y;
				m_z = rhs_.m_z;
				m_w = rhs_.m_w;
			}

			///コピーコンストラクタ
			///[0]i + [1]j + [2]k + [3] 
			TQuaternion(const value_type(&rhs_)[4])
			{
				LOCALPRINTF("TQuaternion array copy constructor\n");
				m_x = rhs_[0];
				m_y = rhs_[1];
				m_z = rhs_[2];
				m_w = rhs_[3];
			}

			///コピーコンストラクタ for expression template
			template<class Lhs, class Op, class Rhs>
			TQuaternion(const TQuaternionExpression<Lhs, Op, Rhs>& obj_)
			:m_x(obj_.x()),
			 m_y(obj_.y()),
			 m_z(obj_.z()),
			 m_w(obj_.w())
			{
				LOCALPRINTF("TQuaternion expression copy constructor\n");
			}

			///コピー
			self_type& operator = (const self_type& rhs_) 
			{
				LOCALPRINTF("TQuaternion copy\n");
				m_x = rhs_.m_x;
				m_y = rhs_.m_y;
				m_z = rhs_.m_z;
				m_w = rhs_.m_w;
				return *this;
			}
				

			///コピー for expression template
			template<class Lhs, class Op, class Rhs>
			self_type& operator = (const TQuaternionExpression<Lhs, Op, Rhs>& rhs_) 
			{
				LOCALPRINTF("TQuaternion expression copy\n");
				
				//A=A*Bのような時に
				//自分の保持している値を直接いじると掛け算の途中で要素の値が入れ替わってしまうため
				//テンポラリオブジェクトに計算結果を保存し最後にスワップ
				self_type obj(rhs_.x(), rhs_.y(), rhs_.z(), rhs_.w());
				*this = obj;
				return *this;
			}

			//-------------------------------------------------------subscript operators	
			value_type& x(void)
			{
				return m_x;
			}			
			const value_type& x(void) const
			{
				return m_x;
			}
			
			value_type& y()
			{
				return m_y;
			}
			const value_type& y() const
			{
				return m_y;
			}

			value_type& z()
			{
				return m_z;
			}
			const value_type& z() const
			{
				return m_z;
			}
			
			value_type& w()
			{
				return m_w;
			}
			const value_type& w() const
			{
				return m_w;
			}


			////-------------------------------------------------------cast operators
			///デバッグ用stingキャスト
			operator std::string () const
			{
				std::stringstream ss;
				ss << "quaternion\n\t";
				ss << "(" << m_x << "*i + " <<m_y << "*j + "<< m_z << "*k) + " << m_w;
				return ss.str();
			}	

			//-------------------------------------------------------unary operators
			///符号反転
			self_type operator - () const
			{
				return self_type(m_y, m_y, m_z, m_w);
			}

			//-------------------------------------------------------compound assignment operators
			///和
			self_type& operator += (const self_type& rhs_) 
			{
				m_x += rhs_.m_x;
				m_y += rhs_.m_y;
				m_z += rhs_.m_z;
				m_w += rhs_.m_w;
				return *this;
			}

			///差
			self_type& operator -= (const self_type& rhs_) 
			{
				m_x -= rhs_.m_x;
				m_y -= rhs_.m_y;
				m_z -= rhs_.m_z;
				m_w -= rhs_.m_w;
				return *this;
			}

			///クォータニオンとの積
			self_type& operator *= (const self_type& rhs_) 
			{
				self_type org(*this);
				m_x = org.m_w*rhs_.m_x + org.m_x*rhs_.m_w + org.m_y*rhs_.m_z - org.m_z*rhs_.m_y;
				m_y = org.m_w*rhs_.m_y + org.m_y*rhs_.m_w - org.m_x*rhs_.m_z + org.m_z*rhs_.m_x;
				m_z = org.m_w*rhs_.m_z + org.m_z*rhs_.m_w + org.m_x*rhs_.m_y - org.m_y*rhs_.m_x;
				m_w = org.m_w*rhs_.m_w - org.m_x*rhs_.m_x - org.m_y*rhs_.m_y - org.m_z*rhs_.m_z;
				return *this;
			}

			///スカラ倍
			self_type& operator *= (const value_type& rhs_)
			{
				m_x *= rhs_;
				m_y *= rhs_;
				m_z *= rhs_;
				m_w *= rhs_;
				return *this;
			}

			//-------------------------------------------------------binary operators
			//クラスで定義する場合、ExpressionTemplateに対応を考えると
			//左手/右手の記述がそれぞれに必要で煩雑になる
			//そのため二項演算子はグローバル関数として定義する
					
	
			//-------------------------------------------------------math operations
			///長さ
			value_type Length() const
			{
				//内部的にはrsqrt(lengthsquared) * lengthsquared
				return function::Sqrt(LengthSquared());
			}
			
			///長さの二乗
			value_type LengthSquared() const
			{
				return m_x*m_x + m_y*m_y + m_z*m_z + m_w*m_w;
			}

			///長さの逆数
			value_type ReciprocalLength() const
			{
				return function::RSqrt(LengthSquared());
			}

			///正規化する
			self_type& Normalize()
			{
				value_type rl = ReciprocalLength();
				
				m_x *= rl;
				m_y *= rl;
				m_z *= rl;
				m_w *= rl;
				return *this;
			}
			
			///共役クォータニオンになる
			///幾何学的には逆方向の回転を表す
			self_type& Conjugate()
			{
				m_x = -m_x;
				m_y = -m_y;
				m_z = -m_z;
				return *this;
			}

			///逆クオータニオンになる
			///元のクオータニオンと積が単位元(1, <0,0,0>)
			///
			///qinv = qconj / ||q||^2
			///
			///回転クオータニオンqrを用いている場合、以下の式で回転できる
			///
			///q' = qr * q * Inv(qr)
			///
			self_type& Inverse()
			{
				value_type ls = ReciprocalLength();
				ls *= ls;
				m_x *= -ls;
				m_y *= -ls;
				m_z *= -ls;
				m_w *= ls;
				return *this;
			}

			///符号逆転する
			///幾何学的には回転方向が変わらず、q = -q である
			void Reverse()
			{
				m_x = -m_x;
				m_y = -m_y;
				m_z = -m_z;
				m_w = -m_w;
			}
			
			///項等クォータニオンになる
			self_type& Identity()
			{
				m_x = 0;
				m_y = 0;
				m_z = 0;
				m_w = 1;
				return *this;
			}
		};
	}
}




#endif