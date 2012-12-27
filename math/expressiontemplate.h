//2010/04/18
//rflab.

#ifndef _HRK_EXPRESSION_TEMPLATE_
#define _HRK_EXPRESSION_TEMPLATE_

#include <cassert>
#include <sstream>

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
		//----------------------------------------------expression functor
		//---------------------------------------------------for vector
		struct CVectorAddition
		{
			template<typename T>
			static T Apply(const T& lhs_, const T& rhs_)
			{
				return lhs_ + rhs_;
			}
		};

		struct CVectorSubtraction
		{
			template<typename T>
			static T Apply(const T& lhs_, const T& rhs_)
			{
				return lhs_ - rhs_;
			}
		};
		
		//---------------------------------------------------for matrix
		struct CMatrixAddition
		{
			template<typename T>
			static T Apply(const T& lhs_, const T& rhs_)
			{
				return lhs_ + rhs_;
			}
		};

		struct CMatrixSubtraction
		{
			template<typename T>
			static typename T Apply(const T& lhs_, const T& rhs_)
			{
				return lhs_ - rhs_;
			}
		};

		//---------------------------------------------------for euler angle
		struct CEulerAnglesAddition
		{
			template<typename T>
			static T Apply(const T& lhs_, const T& rhs_)
			{
				return lhs_ + rhs_;
			}
		};

		struct CEulerAnglesSubtraction
		{
			template<typename T>
			static T Apply(const T& lhs_, const T& rhs_)
			{
				return lhs_ - rhs_;
			}
		};

		//---------------------------------------------------for quaterninon
		struct CQuaternionAddition
		{
			template<typename T>
			static T Apply(const T& lhs_, const T& rhs_)
			{
				return lhs_ + rhs_;
			}
		};

		struct CQuaternionSubtraction
		{
			template<typename T>
			static T Apply(const T& lhs_, const T& rhs_)
			{
				return lhs_ - rhs_;
			}
		};
			
		//----------------------------------------------expression template
		///このクラスは算術演算の構造を保持する<br>
		///算術演算のオーバーロードでオブジェクト自体を返すかわりに<br>
		///このクラスを返し、最後にmatrix、vectorなどで定義した=で式を解決することで<br>
		///式の保存を可能にする。また、一時オブジェクトの生成を最小限に抑えることが出来る。<br>
		///<br>
		///example)<br>
		///rf::math::t_vector3 tmp1(1.0f, 2.0f, 3.0f);<br>
		///rf::math::t_vector3 tmp2(4.0f, 5.0f, 6.0f);<br>
		///rf::math::t_vector3 tmp3(7.0f, 8.0f, 9.0f);<br>
		///rf::math::t_vector3 ret = tmp1 + (tmp2 - tmp3); //この計算でのオーバーヘッド抑えることが出来る<br>
		///cout << ret << endl;<br>　
		///<br>
		template<class Lhs, class Op, class Rhs>
		class TVectorExpression
		{
		public:
			typedef typename Lhs::value_type value_type;
			typedef typename Lhs lhs_type;
			typedef typename Rhs rhs_type;
			static const int dimension = lhs_type::dimension;
			
			typedef TVectorExpression<Lhs, Op, Rhs> self_type;

		public:
			const Lhs& m_l;
			const Rhs& m_r;

		public:
			TVectorExpression(const Lhs& l_, const Rhs& r_):m_l(l_), m_r(r_){}
			
			value_type operator ()(int ix_) const
			{
				return Op::Apply(m_l(ix_), m_r(ix_));
			}

			//-------------------------------------------------------casting
			///デバッグ用stingキャスト
			operator std::string () const
			{
				std::stringstream ss;
				ss << "vector<" << Dim << ">\n\t[";
				for (int i=0;i<Dim;i++)
				{
					ss << std::setw(8) << (*this)(i);
					if(i<Dim-1)
						ss << ", ";
				}
				ss << "]";
				return ss.str();
			}
			//-------------------------------------------------------binary operators
			//クラスで定義する場合、ExpressionTemplateに対応を考えると
			//左手/右手の記述がそれぞれに必要で煩雑になる
			//そのため二項演算子はグローバル関数として定義する
		};


		
		template<class Lhs, class Op, class Rhs>
		class TMatrixExpression
		{
		public:
			typedef typename Lhs::value_type value_type;
			typedef typename Lhs lhs_type;
			typedef typename Rhs rhs_type;
			typedef TMatrixExpression<Lhs, Op, Rhs> self_type;
			static const int row_size    = lhs_type::row_size;
			static const int column_size = lhs_type::column_size;

		public:
			const Lhs& m_l;
			const Rhs& m_r;

		public:
			TMatrixExpression(const Lhs& l_, const Rhs& r_):m_l(l_), m_r(r_){}
			
			value_type operator ()(int row_, int column_) const
			{
				return Op::Apply(m_l(row_, column_), m_r(row_, column_));
			}

			//-------------------------------------------------------casting
			///デバッグ用stingキャスト
			operator std::string ()
			{
				std::stringstream ss;
				ss << "matrix expression<" << row_size << ", "<< column_size << ">\n\t[";
				for (int r=0;r<row_size;r++)
				{
					for (int c=0;c<column_size;c++)
					{
						ss << std::setw(8) << (*this)(r, c); //std::fixed
						if(c<Column-1)
							ss << ", ";
					}
					if(r<Row-1)
						ss << ";\n\t ";
				}
				ss << "]";
				return ss.str();
			}

			//-------------------------------------------------------binary operators
			//クラスで定義する場合、ExpressionTemplateに対応を考えると
			//左手/右手の記述がそれぞれに必要で煩雑になる
			//そのため二項演算子はグローバル関数として定義する
		};
		
		template<class Lhs, class Op, class Rhs>
		class TEulerAnglesExpression
		{
		public:
			typedef typename Lhs::value_type value_type;
			typedef typename Lhs lhs_type;
			typedef typename Rhs rhs_type;
			
			typedef TEulerAnglesExpression<Lhs, Op, Rhs> self_type;

		public:
			const Lhs& m_l;
			const Rhs& m_r;

		public:
			TEulerAnglesExpression(const Lhs& l_, const Rhs& r_):m_l(l_), m_r(r_){}
			
			value_type pitch() const
			{
				return Op::Apply(m_l.x(), m_r.x());
			}
			
			value_type heading() const
			{
				return Op::Apply(m_l.y(), m_r.y());
			}

			value_type bank() const
			{
				return Op::Apply(m_l.z(), m_r.z());
			}

			//-------------------------------------------------------casting
			///デバッグ用stingキャスト
			operator std::string () const
			{
				std::stringstream ss;
				ss << "euler angles expression\n\t";
				ss << "(" << heading() << ", " << pitch() << ","<< bank() << ")";
				return ss.str();
			}

			//-------------------------------------------------------binary operators
			//クラスで定義する場合、ExpressionTemplateに対応を考えると
			//左手/右手の記述がそれぞれに必要で煩雑になる
			//そのため二項演算子はグローバル関数として定義する
		};


		template<class Lhs, class Op, class Rhs>
		class TQuaternionExpression
		{
		public:
			typedef typename Lhs::value_type value_type;
			typedef typename Lhs lhs_type;
			typedef typename Rhs rhs_type;
			
			typedef TQuaternionExpression<Lhs, Op, Rhs> self_type;

		public:
			const Lhs& m_l;
			const Rhs& m_r;

		public:
			TQuaternionExpression(const Lhs& l_, const Rhs& r_):m_l(l_), m_r(r_){}
			
			value_type x() const
			{
				return Op::Apply(m_l.x(), m_r.x());
			}
			
			value_type y() const
			{
				return Op::Apply(m_l.y(), m_r.y());
			}

			value_type z() const
			{
				return Op::Apply(m_l.z(), m_r.z());
			}
			
			value_type w() const
			{
				return Op::Apply(m_l.w(), m_r.w());
			}

			//-------------------------------------------------------casting
			///デバッグ用stingキャスト
			operator std::string () const
			{
				std::stringstream ss;
				ss << "quaternion\n\t";
				ss << "(" << x() << "*i + " << y() << "*j + "<< z() << "*k) + " << w();
				return ss.str();
			}

			//-------------------------------------------------------binary operators
			//クラスで定義する場合、ExpressionTemplateに対応を考えると
			//左手/右手の記述がそれぞれに必要で煩雑になる
			//そのため二項演算子はグローバル関数として定義する
		};

	}
}

#endif