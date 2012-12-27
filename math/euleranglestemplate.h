//2010/04/18
//rflab.

#ifndef _RF_EULER_ANGLES_TEMPLATE_H_
#define _RF_EULER_ANGLES_TEMPLATE_H_

#include <cassert>
#include <sstream>

#include "math.h"
#include "expressiontemplate.h"

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
		///オイラー角
		///pitch→heading→bank
		template<typename T>
		class TEulerAngles
		{
		public:
			T m_heading;
			T m_pitch;
			T m_bank;

		public:
			typedef T value_type;
			typedef TEulerAngles<T> self_type;
		public:

			///デストラクタ
			///virtualはつけない->速度が異様に落ちる
			~TEulerAngles(){}

			///引数なしコンストラクタ
			///初期化も含めて何もしない
			TEulerAngles()
			{
				LOCALPRINTF("TEulerAngles no arg  constructor\n");
			}

			///アングル指定コンストラクタ
			TEulerAngles(T heading_, T pitch_, T bank_)
			{
				LOCALPRINTF("TEulerAngles each member init constructor\n");
				m_heading = heading_;
				m_pitch   = pitch_;
				m_bank    = bank_;
			}


			///コピーコンストラクタ
			TEulerAngles(const self_type& rhs_) 
			{
				LOCALPRINTF("TEulerAngles copy constructor\n");
				m_heading = rhs_.m_heading;
				m_pitch = rhs_.m_pitch;
				m_bank = rhs_.m_bank;
			}

			///コピーコンストラクタ for expression template
			template<class Lhs, class Op, class Rhs>
			TEulerAngles(const TEulerAnglesExpression<Lhs, Op, Rhs>& obj_)
			:m_heading(obj_.heading()),
			 m_pitch(obj_.pitch()),
			 m_bank(obj_.bank())
			{
				LOCALPRINTF("TEulerAngles expression copy constructor\n");
			}

			///コピー
			self_type& operator = (const self_type& rhs_) 
			{
				LOCALPRINTF("TEulerAngles copy\n");
				m_heading = rhs_.m_heading;
				m_pitch   = rhs_.m_pitch;
				m_bank    = rhs_.m_bank;
				return *this;
			}				

			///コピー for expression template
			template<class Lhs, class Op, class Rhs>
			self_type& operator = (const TEulerAnglesExpression<Lhs, Op, Rhs>& rhs_) 
			{
				LOCALPRINTF("TEulerAngles expression copy\n");
				
				//A=A*Bのような時に
				//自分の保持している値を直接いじると掛け算の途中で要素の値が入れ替わってしまうため
				//テンポラリオブジェクトに計算結果を保存し最後にスワップ
				self_type obj(rhs_.heading(), rhs_.pitch(), rhs_.bank());
				*this = obj;
				return *this;
			}

			//-------------------------------------------------------subscript operators	
			value_type& heading(void)
			{
				return m_heading;
			}
				
			const value_type& heading(void) const
			{
				return m_heading;
			}
			
			value_type& pitch()
			{
				return m_pitch;
			}
		
			const value_type& pitch() const
			{
				return m_pitch;
			}
	
			value_type& bank()
			{
				return m_bank;
			}

			const value_type& bank() const
			{
				return m_bank;
			}


			////-------------------------------------------------------cast operators
			///デバッグ用stingキャスト
			operator std::string () const
			{
				std::stringstream ss;
				ss << "euler angles\n\t";
				//ss << m_elem[0] << "+ (" << m_elem[1] << "*i + " << m_elem[2] << "*j + "<< m_elem[3] << "*k)";
				ss << "(" << m_heading << ", " <<m_pitch << ","<< m_bank << ")";
				return ss.str();
			}	

			//-------------------------------------------------------unary operators
			self_type operator - () 
			{
				return self_type(-m_heading, -m_pitch, -m_bank);
			}

			//-------------------------------------------------------compound assignment operators
			///和
			self_type& operator += (const self_type& rhs_) 
			{
				m_heading += rhs_.m_heading;
				m_pitch   += rhs_.m_pitch;
				m_bank    += rhs_.m_bank;
				return *this;
			}

			///差
			self_type& operator -= (const self_type& rhs_) 
			{
				m_heading -= rhs_.m_heading;
				m_pitch   -= rhs_.m_pitch;
				m_bank    -= rhs_.m_bank;
				return *this;
			}

			///スカラ倍
			self_type& operator *= (const value_type& rhs_)
			{
				m_heading *= rhs_;
				m_pitch   *= rhs_;
				m_bank    *= rhs_;
				return *this;
			}
			
			//-------------------------------------------------------binary operators
			//クラスで定義する場合、ExpressionTemplateに対応を考えると
			//左手/右手の記述がそれぞれに必要で煩雑になる
			//そのため二項演算子はグローバル関数として定義する

			//-------------------------------------------------------math operations
			///正準化
			///canonize:（教会が）～を認める
			///いまのところfloat only
			self_type& Canonize()
			{
				//pitchを正準範囲[-π/2, π/2]に変換する
				//[-π, π]の範囲に周回した後に裏表を判定して正準範囲に変換
				math::function::Cycle(&m_pitch, -definition::RF_PI, definition::RF_PI); 
				if (m_pitch < -definition::RF_OVER_2PI)
				{
					m_pitch   -= definition::RF_PI;
					m_heading += definition::RF_PI;
					m_bank    += definition::RF_PI;
				}
				else if(definition::OVER_PI2 < m_pitch)
				{
					m_pitch    = definition::RF_PI - m_pitch;
					m_heading += definition::RF_PI;
					m_bank    += definition::RF_PI;
				}

				//gimbal lockしていなければ
				//headingを正準範囲に周回
				if (std::fabs(m_pitch)>definition::RF_OVER_2PI - definition::RF_EPSILON)
				{
					m_heading += m_bank;
					m_bank = 0.0f;
				}
				else
				{
					math::function::Cycle(&m_bank, -definition::RF_PI, definition::RF_PI); 
				}

				//headingを正準範囲に周回
				math::function::Cycle(&m_heading, -definition::RF_PI, definition::RF_PI); 

				return *this;
			}
			
			///項等オイラー角になる
			self_type& Identity()
			{
				m_heading = 0;
				m_pitch   = 0;
				m_bank    = 0;
				return *this;
			}
		};
	}
}




#endif