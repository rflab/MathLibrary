//2010/04/18
//rflab.

#ifndef _RF_MATRIX_TEMPLATE_
#define _RF_MATRIX_TEMPLATE_

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
		///行列
		template<typename T, int Row, int Column>
		class TMatrix
		{
			//template<class, int, int> friend class TMatrix; ///< for upcast -- TMatrix<任意の型> はTMatrix<T>のfriend
	
		private:
			#if 0 
				廃止 TVector<TVector<T, Column>, Row> m_elem;
			#else
				T m_elem[Row][Column];
			#endif
		public:
			typedef TMatrix<T, Row, Column> self_type;
			typedef T row_type[Row];
			typedef typename T value_type;

			static const int row_size    = Row;
			static const int column_size = Column;

		public:
			
			///デストラクタ
			///virtualはつけない->速度が異様に落ちる
			~TMatrix(){}

			///引数無しコンストラクタ
			///初期化も含めて何もしない
			TMatrix()
			{
				LOCALPRINTF("TMatrix no arg constructor\n");
			}

			
			///1*1要素の場合のコンストラクタ
			TMatrix(const typename T& init_)
			{
				LOCALPRINTF("TMatrix1 initial scalar constructor\n");
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] = init_;
			}

			///2*2要素の場合のコンストラクタ
			TMatrix(
				const T& a00_, const T& a01_,
				const T& a10_, const T& a11_)
			{
				LOCALPRINTF("TMatrix 2*2 args  constructor2\n");
				assert(Row == 2);
				assert(Column == 2);
				m_elem[0][0] = a00_;
				m_elem[0][1] = a01_;
				m_elem[1][0] = a10_;
				m_elem[1][1] = a11_;
			}
		
			///3*3要素の場合のコンストラクタ
			TMatrix(
				const T& a00_, const T& a01_, const T& a02_,
				const T& a10_, const T& a11_, const T& a12_,
				const T& a20_, const T& a21_, const T& a22_)
			{
				LOCALPRINTF("TMatrix 3*3 args constructor3\n");
				assert(Row == 3);
				assert(Column == 3);
				m_elem[0][0] = a00_;
				m_elem[0][1] = a01_;
				m_elem[0][2] = a02_;
				m_elem[1][0] = a10_;
				m_elem[1][1] = a11_;
				m_elem[1][2] = a12_;
				m_elem[2][0] = a20_;
				m_elem[2][1] = a21_;
				m_elem[2][2] = a22_;
			}
		
			///4*4要素の場合のコンストラクタ
			TMatrix(
				const T& a00_, const T& a01_, const T& a02_, const T& a03_,
				const T& a10_, const T& a11_, const T& a12_, const T& a13_,
				const T& a20_, const T& a21_, const T& a22_, const T& a23_,
				const T& a30_, const T& a31_, const T& a32_, const T& a33_)
			{
				LOCALPRINTF("TMatrix 4*4 args constructor4\n");
				assert(Row == 4);
				assert(Column == 4);
				m_elem[0][0] = a00_;
				m_elem[0][1] = a01_;
				m_elem[0][2] = a02_;
				m_elem[0][3] = a03_;
				m_elem[1][0] = a10_;
				m_elem[1][1] = a11_;
				m_elem[1][2] = a12_;
				m_elem[1][3] = a13_;
				m_elem[2][0] = a20_;
				m_elem[2][1] = a21_;
				m_elem[2][2] = a22_;
				m_elem[2][3] = a23_;
				m_elem[3][0] = a30_;
				m_elem[3][1] = a31_;
				m_elem[3][2] = a32_;
				m_elem[3][3] = a33_;
			}

			/////      	コンストラクタテンプレート
			/////D3DXMATRIXからのコピーのため
			/////危険なことかも知れませんがよく分かりません
			//template<typename S>
			//TMatrix(const S& rhs_) 
			//{
			//}	
				
			/////      	配列によるコピーコンストラクタ
			//TMatrix(const T (&init_)[Row][Column])
			//{
			//	LOCALPRINTF("TMatrix array copy constructor\n");
			//	for(int r=0; r<Row; r++)
			//	for(int c=0; c<Column; c++)
			//		m_elem[r][c] = init_[r][c];
			//}
			
			///配列ポインタコピーコンストラクタ
			///代入演算子も用意しないとコピーコンストラクタで一時オブジェクト作った後、
			///代入が行われる？
			TMatrix(const value_type* rhs_) 
			{
				LOCALPRINTF("TMatrix array pointer copy constructor\n");

				int cnt = 0;
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] = rhs_[cnt++];
			}	
			
			///コピーコンストラクタ
			TMatrix(const self_type& rhs_) 
			{
				LOCALPRINTF("TMatrix copy constructor\n");
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] = rhs_(r, c);
			}

			//グローバル関数として提供することにした
			/////      	コピーコンストラクタ
			/////クォータニオンから同等の回転行列を作成
			//template<template<typename>class Quaternion>
			//TMatrix(const Quaternion<value_type>& rhs_) 
			//{
			//	LOCALPRINTF("TMatrix copy constructor quaternion\n");
			//	assert(Row == 4);
			//	assert(Column == 4);
			//	Quaternion<value_type> q2(rhs_);
			//	q2 += q2;
			//	Quaternion<value_type> qw(rhs_.x()*rhs_.x(), rhs_.y()*rhs_.y(), rhs_.z()*rhs_.z(), rhs_.w()*rhs_.w()); //２乗
			//	qw += qw;
			//	value_type _2xy = q2.x() * rhs_.y();
			//	value_type _2yz = q2.y() * rhs_.z();
			//	value_type _2zx = q2.x() * rhs_.z();
			//	value_type _2xw = rhs_.x() * q2.w();
			//	value_type _2yw = rhs_.y() * q2.w();
			//	value_type _2zw = rhs_.z() * q2.w();
			//	
			//	m_elem[0][0] = 1 - qw.y() - qw.z();
			//	m_elem[0][1] = _2xy + _2zw;
			//	m_elem[0][2] = _2zx - _2yw;
			//	m_elem[0][3] = 0;
			//	m_elem[1][0] = _2xy - _2zw;;
			//	m_elem[1][1] = 1 - qw.z() - qw.x();
			//	m_elem[1][2] = _2yz + _2xw;
			//	m_elem[1][3] = 0;
			//	m_elem[2][0] = _2zx + _2yw;
			//	m_elem[2][1] = _2yz - _2xw;
			//	m_elem[2][2] = 1 - qw.x() - qw.y();
			//	m_elem[2][3] = 0;
			//	m_elem[3][0] = 0;
			//	m_elem[3][1] = 0;
			//	m_elem[3][2] = 0;
			//	m_elem[3][3] = 1;
			//}

			///コピーコンストラクタ for expression template
			template<class Lhs, class Op, class Rhs>
			TMatrix(const TMatrixExpression<Lhs, Op, Rhs>& obj_)
			{
				LOCALPRINTF("TVector expression evaluation copy constructor\n");
				
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
				{
					m_elem[r][c] = obj_(r, c);
				}
			}

			///コピー
			self_type& operator = (const self_type& rhs_) 
			{
				LOCALPRINTF("TMatrix copy\n");
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] = rhs_(r, c);
				return *this;
			}

			///配列ポインタコピー
			self_type& operator = (const value_type* rhs_) 
			{
				LOCALPRINTF("TMatrix array pointer copy\n");
				int cnt = 0;
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] = rhs_[cnt++];
				return *this;
			}
				
			///コピー
			template<class Lhs, class Op, class Rhs>
			TMatrix& operator = (const TMatrixExpression<Lhs, Op, Rhs>& rhs_) 
			{
				LOCALPRINTF("TVector expression evaluation copy\n");

				//乗算はexpressionの中で解決することにしたので
				//一時オブジェクトの生成は必要なくなった
				#if 0
					//自分の保持している値を直接いじると掛け算とかで何回も同じ値を参照する
					//計算する値と結果を保存する値が同じ場合に
					//なので、一時オブジェクトに書き込んでスワップする
					self_type tmp;
					for(int r=0; r<Row; r++)
					for(int c=0; c<Column; c++)
					{
						tmp.m_elem[r][c] = rhs_(r, c);
					}
					*this = tmp;
				#else
					for(int r=0; r<Row; r++)
					for(int c=0; c<Column; c++)
					{
						m_elem[r][c] = rhs_(r, c);
					}
				#endif
				return *this;
			}

			//-------------------------------------------------------casting
			///配列化キャスト
			//operator LPTSTR ();のようなもの
			operator T* ()
			{
				return (T*)m_elem;
			}

			///配列化キャスト
			operator const T* () const
			{
				return m_elem;
			}	
			
			///デバッグ用stingキャスト
			operator std::string ()
			{
				std::stringstream ss;
				ss << "matrix<" << Row << ", "<< Column << ">\n\t[";
				for (int r=0;r<Row;r++)
				{
					for (int c=0;c<Column;c++)
					{
						ss << std::setw(8) << m_elem[r][c]; //std::fixed
						if(c<Column-1)
							ss << ", ";
					}
					if(r<Row-1)
						ss << ";\n\t ";
				}
				ss << "]";
				return ss.str();
			}

			//-------------------------------------------------------subscript operators
			///()による配列アクセス
			value_type& operator()(int row_, int colmun_)
			{
				return m_elem[row_][colmun_];
			}
			const value_type& operator()(int row_, int colmun_) const
			{
				return m_elem[row_][colmun_];
			}

			//-------------------------------------------------------unary operators
			self_type operator - () const
			{
				self_type result;
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					result(r, c) = -m_elem[r][c];
				return result;
			}

			//-------------------------------------------------------compound assignment operators
			///和
			self_type& operator += (const self_type& rhs_) 
			{
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] += rhs_.m_elem[r][c];
				return *this;
			}

			///差
			self_type& operator -= (const self_type& rhs_) 
			{
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] -= rhs_.m_elem[r][c];
				return *this;
			}
			
			///積
			///ただし、行列はColumn×Row（要素数反転）の正方行列に限る
			template <template<typename S, int SRow, int SColumn>class Rhs>
			self_type& operator *= (const Rhs<T, Column, Row>& rhs_) 
			{
				self_type org(*this);

				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
				{
					m_elem[r][c] = org.m_elem[r][0] * rhs_.m_elem[0][r];
					for (int i=1; i<Column ;i++)
					{
						m_elem[r][c] += org.m_elem[r][i] * rhs_.m_elem[i][r];
					}
				}
				return *this;
			}

			//-------------------------------------------------------binary operators
			//クラスで定義する場合、ExpressionTemplateに対応を考えると
			//左手/右手の記述がそれぞれに必要で煩雑になる
			//そのため二項演算子はグローバル関数として定義する
					
	
			//-------------------------------------------------------math operations
			///転置
			void Transpose()
			{
				assert(Column == Row);
				for(int c=1; c<Column; c++)
				for(int r=0; r<c; r++)
					std::swap(m_elem[r][c], m_elem[c][r]);
			}

			///単位行列になる
			void Identity()
			{
				const value_type unit(1);
				const value_type zero(0);
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] = r == c ? unit : zero;
			}
			
			///零行列になる
			void Zero()
			{
				const value_type zero(0);
				for(int r=0; r<Row; r++)
				for(int c=0; c<Column; c++)
					m_elem[r][c] = zero;
			}
		};	
	}
}

#endif