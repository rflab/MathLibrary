//2010/04/18
//rflab.

#ifndef _RF_VECTOR_TEMPLATE_H_
#define _RF_VECTOR_TEMPLATE_H_

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
		///ベクトル
		///directxにあわせて行ベクトル（横向き）
		template<typename T, int Dim>
		class TVector
		{
		private:
			T m_elem[Dim];
			
		public:
			typedef T value_type;
			typedef TVector<T, Dim> self_type;
			static const int dimension = Dim;
		public:
			
			///デストラクタ
			///virtualはつけない->速度が異様に落ちる
			~TVector(){}

			///引数無しコンストラクタ
			///初期化も含めて何もしない
			TVector()
			{
				LOCALPRINTF("TVector no arg constructor\n");
			}
						
			#if 0
				///1要素の場合のコンストラクタ、スカラーで初期化
				explicit TVector(const T& init_)
				{
					LOCALPRINTF("TVector scalar constructor\n");
					for(int i=0; i<Dim; i++)
						m_elem[i] = init_;
				}
			#endif
			
			///2要素の場合のコンストラクタ
			TVector(const T& a0_, const T& a1_)
			{
				LOCALPRINTF("TVector 2 args constructor\n");
				assert(Dim == 2);
				m_elem[0] = a0_;
				m_elem[1] = a1_;
			}
		
			///3要素の場合のコンストラクタ
			TVector(const T& a0_, const T& a1_, const T& a2_)
			{
				LOCALPRINTF("TVector 3 args constructor\n");
				assert(Dim == 3);
				m_elem[0] = a0_;
				m_elem[1] = a1_;
				m_elem[2] = a2_;
			}

			///4要素の場合のコンストラクタ
			TVector(const T& a0_, const T& a1_, const T& a2_, const T& a3_)
			{
				LOCALPRINTF("TVector 4 args constructor\n");
				assert(Dim == 4);
				m_elem[0] = a0_;
				m_elem[1] = a1_;
				m_elem[2] = a2_;
				m_elem[3] = a3_;
			}
			
			///5要素の場合のコンストラクタ
			TVector(const T& a0_, const T& a1_, const T& a2_, const T& a3_, const T& a4_)
			{
				LOCALPRINTF("TVector 5 args constructor\n");
				assert(Dim == 5);
				m_elem[0] = a0_;
				m_elem[1] = a1_;
				m_elem[2] = a2_;
				m_elem[3] = a3_;
				m_elem[4] = a4_;
			}
			
			///6要素の場合のコンストラクタ
			TVector(const T& a0_, const T& a1_, const T& a2_, const T& a3_, const T& a4_,  const T& a5_)
			{
				LOCALPRINTF("TVector 6 args constructor\n");
				assert(Dim == 6);
				m_elem[0] = a0_;
				m_elem[1] = a1_;
				m_elem[2] = a2_;
				m_elem[3] = a3_;
				m_elem[4] = a4_;
				m_elem[5] = a5_;
			}

			///7要素の場合のコンストラクタ
			TVector(const T& a0_, const T& a1_, const T& a2_, const T& a3_, const T& a4_, const T& a5_, const T& a6_)
			{
				LOCALPRINTF("TVector 7 args constructor\n");
				assert(Dim == 6);
				m_elem[0] = a0_;
				m_elem[1] = a1_;
				m_elem[2] = a2_;
				m_elem[3] = a3_;
				m_elem[4] = a4_;
				m_elem[5] = a5_;
				m_elem[6] = a6_;
			}
			
				
			#if 0
				D3DXVECTOR3などからの変換用に作ったが危険なにおいがしたので廃止
				/////      	コンストラクタテンプレート(変換の定義)
				//template<typename S>
				//TVector(const S& rhs_) 
				//{
				//	LOCALPRINTF("TVector constructor template\n");
				//	for(int i=0; i<Dim; i++)
				//		m_elem[i] = rhs_[i];
				//}
			#endif

			///配列ポインタコピーコンストラクタ
			///代入演算子も用意しないとコピーコンストラクタで一時オブジェクト作った後、
			///代入が行われる？
			#if 1
				TVector(const value_type* rhs_) 
				{
					LOCALPRINTF("TVector array pointer copy constructor\n");
		
					for(int i=0; i<Dim; i++)
						m_elem[i] = rhs_[i];
				}
			#else
				TVector(const T (&init_)[Dim])
				{
					LOCALPRINTF("TVector array copy constructor\n");
					for(int i=0; i<Dim; i++)
						m_elem[i] = init_[i];
				}
			#endif

			///コピーコンストラクタ
			TVector(const self_type& rhs_) 
			{
				LOCALPRINTF("TVector copy constructor\n");
				for(int i=0; i<Dim; i++)
					m_elem[i] = rhs_[i];
			}

			///コピーコンストラクタ for expression template
			template<class Lhs, class Op, class Rhs>
			TVector(const TVectorExpression<Lhs, Op, Rhs>& obj_)
			{
				LOCALPRINTF("TVector expression evaluation copy constructor\n");
				
				for(int i=0; i<Dim; i++)
				{
					m_elem[i] = obj_(i);
				}
			}

			///コピー
			self_type& operator = (const self_type& rhs_) 
			{
				LOCALPRINTF("TVector copy\n");
				for(int i=0; i<Dim; i++)
					m_elem[i] = rhs_.m_elem[i];
				return *this;
			}

			#if 0
				///配列コピー
				self_type& operator = (const T (&init_)[Dim]) 
				{
					LOCALPRINTF("TVector array copy\n");
					for(int i=0; i<Dim; i++)
						m_elem[i] = init_[i];
					return *this;
				}
			#else
				///配列ポインタコピー
				self_type& operator = (const value_type* rhs_) 
				{
					LOCALPRINTF("TVector array pointer copy\n");
					for(int i=0; i<Dim; i++)
						m_elem[i] = rhs_[i];
					return *this;
				}
			#endif

			///コピー for expression template
			///基本的にこれよりもコピーコンストラクタをつかったほうが得
			template<class Lhs, class Op, class Rhs>
			self_type& operator = (const TVectorExpression<Lhs, Op, Rhs>& rhs_) 
			{
				LOCALPRINTF("TVector expression evaluation copy\n");
				
				//乗算はexpressionの中で解決することにしたので
				//一時オブジェクトの生成は必要なくなった
				#if 0
					//A=A*Bのような時に
					//自分の保持している値を直接いじると掛け算の途中で要素の値が入れ替わってしまうため
					//テンポラリオブジェクトに計算結果を保存し最後にスワップ
					self_type obj();
					for(int i=0; i<Dim; i++)
					{
						obj.m_elem[i] = rhs_(i);
					}
					*this = obj;
				#else
					for(int i=0; i<Dim; i++)
					{
						m_elem[i] = rhs_(i);
					}
				#endif

				return *this;
			}

			//-------------------------------------------------------casting
			///配列化キャスト
			operator T* ()
			{
				return m_elem;
			}

			///配列化キャスト
			operator const T* () const
			{
				return m_elem;
			}

			///デバッグ用stingキャスト
			operator std::string () const
			{
				std::stringstream ss;
				ss << "vector<" << Dim << ">\n\t[";
				for (int i=0;i<Dim;i++)
				{
					ss << std::setw(8) << m_elem[i];
					if(i<Dim-1)
						ss << ", ";
				}
				ss << "]";
				return ss.str();
			}


			//-------------------------------------------------------subscript operators			
			///()による配列アクセス
			value_type& operator()(int ix_)
			{
				return m_elem[ix_];
			}
			const value_type& operator()(int ix_) const
			{
				return m_elem[ix_];
			}

			///[]による配列アクセス
			value_type& operator[](int ix_)
			{
				return m_elem[ix_];
			}
			const value_type& operator[](int ix_) const
			{
				return m_elem[ix_];
			}
			
			//-------------------------------------------------------unary operators
			self_type operator - () const
			{
				self_type result;
				for (int i = 0; i < Dim; ++i)
					result(i) = -m_elem[i];
				return result;
			}

			//-------------------------------------------------------compound assignment operators
			///和
			self_type& operator += (const self_type& rhs_) 
			{
				for (int i = 0; i < Dim; ++i)
					m_elem[i] += rhs_.m_elem[i];
				return *this;
			}

			
			///差
			self_type& operator -= (const self_type& rhs_) 
			{
				for (int i = 0; i < Dim; ++i)
					m_elem[i] -= rhs_.m_elem[i];
				return *this;
			}

			///行列との積
			///ただし、行列はDim×Dimの正方行列に限る
			template <template<typename S, int Row, int Column>class Matrix>
			self_type& operator *= (const Matrix<T, Dim, Dim>& rhs_) 
			{
				self_type org(*this);

				for (int i = 0; i < Dim; ++i)
				{
					m_elem[i] = org(0) * rhs_(0, i);
					for (int r = 1; r < Dim; ++r)
						m_elem[i] += org(r) * rhs_(r, i);
				}
				return *this;
			}

			///スカラーとの積
			self_type& operator *= (const value_type& rhs_)
			{
				for (int i = 0; i < Dim; ++i)
					m_elem[i] *= rhs_;
				return *this;
			}
			
			//-------------------------------------------------------math operators
			///長さ（長さの逆数×長さの二乗)
			value_type Length() const
			{
				//内部的にはrsqrt(lengthsquared) * lengthsquared
				return function::Sqrt(LengthSquared());
			}

			///長さの逆数
			value_type ReciprocalLength() const
			{
				return function::RSqrt(LengthSquared());
			}

			///長さの二乗
			value_type LengthSquared() const
			{
				value_type tmp;
				tmp = (m_elem[0] * m_elem[0]);
				for (int i = 1; i < Dim; ++i)
					tmp += (m_elem[i] * m_elem[i]);
				return tmp;
			}

			///正規化する
			void Normalize()
			{
				T r = ReciprocalLength();
				
				for (int i = 0; i < Dim; ++i)
					m_elem[i] *= r;
			}

			///符号反転する
			void Reverse()
			{
				for (int i = 0; i < Dim; ++i)
					m_elem[i] = -m_elem[i];
			}
		};
	}
}

#endif