//2010/04/18
//rflab.

#ifndef _RF_EXPRESSION_OPERATORS_H_
#define _RF_EXPRESSION_OPERATORS_H_

#include "expressiontemplate.h"
#include "vectortemplate.h"
#include "matrixtemplate.h"
#include "quaterniontemplate.h"
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
		///定義されている演算は以下の通り 
		///-------------------------------vector
		///vector + expression 
		///expression + vector
		///vector + vector
		///expression + expression
		///vector - expression
		///expression - vector
		///vector - vector
		///expression - expression
		///vector * scalar
		///scalar * vector
		///-------------------------------matrix
		///matrix + expression
		///expression + matrix
		///matrix + matrix
		///expression + expression
		///matrix - expression
		///expression - matrix
		///matrix - matrix
		///expression - expression
		///matrix * expression
		///expression * matrix
		///matrix * matrix
		///expression * expression
		///matrix * scalar
		///scalar * matrix
		///expression * scalar
		///scalar * expression
		///-------------------------------quaternion
		///quaternion + expression
		///expression + quaternion
		///quaternion + quaternion
		///expression + expression
		///quaternion - expression
		///expression - quaternion
		///quaternion - quaternion
		///expression - expression
		///quaternion * expression
		///expression * quaternion
		///quaternion * quaternion
		///expression * expression
		///quaternion * scalar
		///scalar * quaternion
		///-------------------------------combinations
		///vector * matrix
		///vector * matrix_expression				//unsupported
		///vector_expression * matrix				//unsupported
		///vector_expression * matrix_expression	//unsupported

		//------------------------------------------------compound assignment operators for Matrix
		///vector + expression 
		template<typename Type, int Dim, class Rhs>
		inline TVectorExpression<TVector<Type, Dim>, CVectorAddition, Rhs>
			operator + (const TVector<Type, Dim>& lhs_, const Rhs& rhs_)
		{
			return TVectorExpression<TVector<Type, Dim>, CVectorAddition, Rhs>(lhs_, rhs_);
		}

		///expression + vector
		template<typename Type, int Dim, class Lhs>
		inline TVectorExpression<Lhs, CVectorAddition, TVector<Type, Dim> >
			operator + (const Lhs& lhs_, const TVector<Type, Dim>& rhs_)
		{
			return TVectorExpression<Lhs, CVectorAddition, TVector<Type, Dim> >(lhs_, rhs_);
		}

		///vector + vector
		template<typename Type, int Dim>
		inline TVectorExpression<TVector<Type, Dim>, CVectorAddition, TVector<Type, Dim> >
			operator + (const TVector<Type, Dim>& lhs_, const TVector<Type, Dim>& rhs_)
		{
			return TVectorExpression<TVector<Type, Dim>, CVectorAddition, TVector<Type, Dim> >(lhs_, rhs_);
		}

		///expression + expression
		template<class LLhs, typename LExp, class LRhs, class RLhs, typename RExp, class RRhs>
		inline TMatrixExpression<
			TVectorExpression<LLhs, LExp, LRhs>,
			CMatrixAddition,
			TVectorExpression<RLhs, RExp, RRhs> >
		operator + (
			const TVectorExpression<LLhs, LExp, LRhs>& lhs_,
			const TVectorExpression<RLhs, RExp, RRhs>& rhs_)
		{
			return TMatrixExpression<
				TVectorExpression<LLhs, LExp, LRhs>,
				CMatrixAddition,
				TVectorExpression<RLhs, RExp, RRhs> >(lhs_, rhs_);
		}

		///vector - expression
		template<typename Type, int Dim, class Rhs>
		inline TVectorExpression<TVector<Type, Dim>, CVectorSubtraction, Rhs>
			operator - (const TVector<Type, Dim>& lhs_, const Rhs& rhs_)
		{
			return TVectorExpression<TVector<Type, Dim>, CVectorSubtraction, Rhs>(lhs_, rhs_);
		}

		///expression - vector
		template<typename Type, int Dim, class Lhs>
		inline TVectorExpression<Lhs, CVectorSubtraction, TVector<Type, Dim> >
			operator - (const Lhs& lhs_, const TVector<Type, Dim>& rhs_)
		{
			return TVectorExpression<Lhs, CVectorSubtraction, TVector<Type, Dim> >(lhs_, rhs_);
		}

		///vector - vector
		template<typename Type, int Dim>
		inline TVectorExpression<TVector<Type, Dim>, CVectorSubtraction, TVector<Type, Dim> >
			operator - (const TVector<Type, Dim>& lhs_, const TVector<Type, Dim>& rhs_)
		{
			return TVectorExpression<TVector<Type, Dim>, CVectorSubtraction, TVector<Type, Dim> >(lhs_, rhs_);
		}

		///expression - expression
		template<class LLhs, typename LExp, class LRhs, class RLhs, typename RExp, class RRhs>
		inline TMatrixExpression<
			TVectorExpression<LLhs, LExp, LRhs>,
			CMatrixSubtraction,
			TVectorExpression<RLhs, RExp, RRhs> >
		operator - (
			const TVectorExpression<LLhs, LExp, LRhs>& lhs_,
			const TVectorExpression<RLhs, RExp, RRhs>& rhs_)
		{
			return TMatrixExpression<
				TVectorExpression<LLhs, LExp, LRhs>,
				CMatrixSubtraction,
				TVectorExpression<RLhs, RExp, RRhs> >(lhs_, rhs_);
		}

		///vector * scalar
		///乗算の場合は計算結果を出しておかないと掛け算の繰り返しで損をするのでオブジェクト返し
		template<typename Type, int Dim>
		inline TVector<Type, Dim>
			operator * (
				const TVector<Type, Dim>& lhs_,
				const typename TVector<Type, Dim>::value_type& rhs_)
		{	
			TVector<Type, Dim> obj(lhs_);
			for (int i=0; i<Dim; i++)
			{
				obj(i) = lhs_(i) * rhs_;
			}
			return obj;
		}

		///scalar * vector
		///乗算の場合は計算結果を出しておかないと掛け算の繰り返しで損をするのでオブジェクト返し
		template<typename Type, int Dim>
		inline TVector<Type, Dim>
			operator * (
				const typename TVector<Type, Dim>::value_type& lhs_,
				const TVector<Type, Dim>& rhs_)
		{	
			TVector<Type, Dim> obj;
			for (int i=0; i<Dim; i++)
			{
				obj(i) = lhs_ * rhs_(i);
			}
			return obj;
		}

		//------------------------------------------------compound assignment operators for Matrix
		///matrix + expression
		template<typename Type, int Row, int Colmun, class Rhs>
		inline TMatrixExpression<TMatrix<Type, Row, Colmun>, CMatrixAddition, Rhs>
			operator + (const TMatrix<Type, Row, Colmun>& lhs_, const Rhs& rhs_)
		{
			return TMatrixExpression<TMatrix<Type, Row, Colmun>, CMatrixAddition, Rhs>(lhs_, rhs_);
		}

		///expression + matrix
		template<typename Type, int Row, int Colmun, class Lhs>
		inline TMatrixExpression<Lhs, CMatrixAddition, TMatrix<Type, Row, Colmun> >
		operator + (const Lhs& lhs_, const TMatrix<Type, Row, Colmun>&rhs_)
		{
			return TMatrixExpression<Lhs, CMatrixAddition, TMatrix<Type, Row, Colmun> >(lhs_, rhs_);
		}

		///matrix + matrix
		template<typename Type, int Row, int Colmun>
		inline TMatrixExpression<
			TMatrix<Type, Row, Colmun>,
			CMatrixAddition,
			TMatrix<Type, Row, Colmun> >
		operator + (
			const TMatrix<Type, Row, Colmun>& lhs_,
			const TMatrix<Type, Row, Colmun>& rhs_)
		{
			return TMatrixExpression<
				TMatrix<Type, Row, Colmun>,
				CMatrixAddition,
				TMatrix<Type, Row, Colmun> >(lhs_, rhs_);
		}
		
		///expression + expression
		template<class LLhs, typename LExp, class LRhs, class RLhs, typename RExp, class RRhs>
		inline TMatrixExpression<
			TMatrixExpression<LLhs, LExp, LRhs>,
			CMatrixAddition,
			TMatrixExpression<RLhs, RExp, RRhs> >
		operator + (
			const TMatrixExpression<LLhs, LExp, LRhs>& lhs_,
			const TMatrixExpression<RLhs, RExp, RRhs>& rhs_)
		{
			return TMatrixExpression<
				TMatrixExpression<LLhs, LExp, LRhs>,
				CMatrixAddition,
				TMatrixExpression<RLhs, RExp, RRhs> >(lhs_, rhs_);
		}

		///matrix - expression
		template<typename Type, int Row, int Colmun, class Rhs>
		inline TMatrixExpression<TMatrix<Type, Row, Colmun>, CMatrixSubtraction, Rhs>
		operator - (const TMatrix<Type, Row, Colmun>&, const Rhs& rhs_)
		{
			return TMatrixExpression<TMatrix<Type, Row, Colmun>, CMatrixSubtraction, Rhs>(lhs_, rhs_);
		}

		///expression - matrix
		template<typename Type, int Row, int Colmun, class Lhs>
		inline TMatrixExpression<Lhs, CMatrixSubtraction, TMatrix<Type, Row, Colmun> >
		operator - (const Lhs& lhs_, const TMatrix<Type, Row, Colmun>& rhs_)
		{
			return TMatrixExpression<Lhs, CMatrixSubtraction, TMatrix<Type, Row, Colmun> >(lhs_, rhs_);
		}

		///matrix - matrix
		template<typename Type, int Row, int Colmun>
		inline TMatrixExpression<
			TMatrix<Type, Row, Colmun>,
			CMatrixSubtraction,
			TMatrix<Type, Row, Colmun> >
		operator - (
			const TMatrix<Type, Row, Colmun>& lhs_,
			const TMatrix<Type, Row, Colmun>& rhs_)
		{
			return TMatrixExpression<
				TMatrix<Type, Row, Colmun>, 
				CMatrixSubtraction,
				TMatrix<Type, Row, Colmun> >(lhs_, rhs_);
		}

		///expression - expression
		template<class LLhs, typename LExp, class LRhs, class RLhs, typename RExp, class RRhs>
		inline TMatrixExpression<
			TMatrixExpression<LLhs, LExp, LRhs>,
			CMatrixSubtraction,
			TMatrixExpression<RLhs, RExp, RRhs> >
		operator - (
			const TMatrixExpression<LLhs, LExp, LRhs>& lhs_,
			const TMatrixExpression<RLhs, RExp, RRhs>& rhs_)
		{
			return TMatrixExpression<
				TMatrixExpression<LLhs, LExp, LRhs>,
				CMatrixSubtraction,
				TMatrixExpression<RLhs, RExp, RRhs> >(lhs_, rhs_);
		}

		///matrix * expression
		///*=でOKかしら?
		///乗算の場合は計算結果を出しておかないと掛け算の繰り返しで損をするのでオブジェクト返し
		template<typename Type, int LRow, int LColmun, typename Rhs>
		inline TMatrix<Type, LRow, Rhs::column_size>
			operator * (
				const TMatrix<Type, LRow, LColmun>& lhs_,
				const Rhs& rhs_) 
		{
			//evaluate expression
			TMatrix<Type, Rhs::row_size, Rhs::column_size> rhs(rhs_);

			//create temporary object for reduce matrix multiplication overhead
			TMatrix<Type, LRow, Rhs::column_size> tmp;

			for(int r=0; r<LRow; r++) 
			for(int c=0; c<Rhs::column_size; c++)
			{
				tmp(r, c) = lhs_(r, 0) * rhs(0, c);
				for (int i=1; i<LColmun ;i++)
				{
					tmp(r, c) += lhs_(r, i) * rhs(i, c);
				}
			}
			return tmp;
		}

		///expression * matrix
		///*=でOKかしら?
		///乗算の場合は計算結果を出しておかないと掛け算の繰り返しで損をするのでオブジェクト返し
		template<typename Lhs, typename Type, int RRow, int RColmun>
		inline TMatrix<Type, Lhs::row_size, RColmun>
			operator * (
				const Lhs& lhs_,
				const TMatrix<Type, RRow, RColmun>& rhs_) 
		{
			//evaluate expression
			TMatrix<Type, Lhs::row_size, Lhs::column_size> lhs(lhs_);

			//create temporary object for reduce matrix multiplication overhead
			TMatrix<Type, Lhs::row_size, RColmun> tmp;

			for(int r=0; r<Lhs::row_size; r++) 
			for(int c=0; c<RColmun; c++)
			{
				tmp(r, c) = lhs(r, 0) * rhs_(0, c);
				for (int i=1; i<Lhs::column_size ;i++)
				{
					tmp(r, c) += lhs(r, i) * rhs_(i, c);
				}
			}
			return tmp;
		}

		///matrix * matrix
		///*=でOKかしら?
		///乗算の場合は計算結果を出しておかないと掛け算の繰り返しで損をするのでオブジェクト返し
		template<typename Type, int LRow, int CommonDim, int RColmun>
		inline TMatrix<Type, LRow, RColmun>
			operator * (
				const TMatrix<Type, LRow, CommonDim>& lhs_,
				const TMatrix<Type, CommonDim, RColmun>& rhs_)
		{
			//create temporary object for reduce matrix multiplication overhead
			TMatrix<Type, LRow, RColmun> tmp;

			for(int r=0; r<LRow; r++) 
			for(int c=0; c<RColmun; c++)
			{
				tmp(r, c) = lhs_(r, 0) * rhs_(0, c);
				for (int i=1; i<CommonDim ;i++)
				{
					tmp(r, c) += lhs_(r, i) * rhs_(i, c);
				}
			}
			return tmp;
		}

		///expression * expression
		///*=でOKかしら?
		///乗算の場合は計算結果を出しておかないと掛け算の繰り返しで損をするのでオブジェクト返し
		template<class LLhs, typename LExp, class LRhs, class RLhs, typename RExp, class RRhs>
		inline TMatrix<typename LLhs::value_type, LLhs::row_size, RLhs::column_size>
			operator * (
				const TMatrixExpression<LLhs, LExp, LRhs>& lhs_,
				const TMatrixExpression<RLhs, RExp, RRhs>& rhs_)
		{
			//create temporary object for reduce matrix multiplication overhead
			TMatrix<LLhs::value_type, LRow, LColmun> lhs(lhs_);
			TMatrix<LLhs::value_type, RRow, RColmun> rhs(rhs_);

			//create temporary object for reduce matrix multiplication overhead
			TMatrix<LLhs::value_type, LLhs::row_size, RLhs::column_size> tmp;

			for(int r=0; r<LLhs::row_size; r++) 
			for(int c=0; c<RLhs::column_size; c++)
			{
				tmp(r, c) = lhs(r, 0) * rhs(0, c);
				for (int i=1; i<LLhs::column_size;i++)
				{
					tmp(r, c) += lhs(r, i) * rhs(i, c);
				}
			}
			return tmp;
		}
		///matrix * scalar
		///*=でOKかしら?
		///これはexpression templateに出来る気がするが面倒なので暫定で一時オブジェクトを作成する
		template<typename Type, int Row, int Column>
		inline TMatrix<Type, Row, Column>
			operator * (
				const TMatrix<Type, Row, Column>& lhs_,
				const typename TMatrix<Type, Row, Column>::value_type&  rhs_)
		{	
			TMatrix<Type, Row, Column> obj;

			for(int r=0; r<Row; r++)
			for(int c=0; c<Column; c++)
			{
				obj(r, c) = lhs_(r, c) * rhs_;
			}
			return obj;
		}

		///scalar * matrix
		///*=でOKかしら?
		///これはexpression templateに出来る気がするが面倒なので暫定で一時オブジェクトを作成する
		template<typename Type, int Row, int Column>
		inline TMatrix<Type, Row, Column>
			operator * (
				const typename TMatrix<Type, Row, Column>::value_type& lhs_,
				const TMatrix<Type, Row, Column>& rhs_)
		{	
			TMatrix<Type, Row, Column> obj;

			for(int r=0; r<Row; r++)
			for(int c=0; c<Column; c++)
			{
				obj(r, c) = lhs_ * rhs_(r, c);
			}
			return obj;
		}

		
		///expression * scalar
		///面倒なのでオブジェクト返し、いつかExpression Templateにする
		template<class LLhs, typename LExp, class LRhs>
		inline TMatrix<typename LLhs::value_type, LLhs::row_size, LLhs::column_size>
			operator * (
				const TMatrixExpression<LLhs, LExp, LRhs>& lhs_,
				const typename LLhs::value_type& rhs_)
		{	
			TMatrix<LLhs::value_type, LLhs::row_size, LRhs::column_size> result(lhs_);
			
			for(int r=0; r<LLhs::row_size; r++) 
			for(int c=0; c<LLhs::column_size; c++)
			{
				result(r, c) = lhs_(r, c) * rhs_;
			}
			return result;
		}	
		///scalar * expression
		///面倒なのでオブジェクト返し、いつかExpression Templateにする
		template<class RLhs, typename RExp, class RRhs>
		inline TMatrix<typename RLhs::value_type, RLhs::row_size, RRhs::column_size>
			operator * (
				const typename RLhs::value_type& rhs_,
				const TMatrixExpression<RLhs, RExp, RRhs>& lhs_)
		{	
			TMatrix<RLhs::value_type, RLhs::row_size, RLhs::column_size> result(lhs_);
			
			for(int r=0; r<RLhs::row_size; r++) 
			for(int c=0; c<RLhs::column_size; c++)
			{
				result(r, c) = lhs_(r, c) * rhs_;
			}
			return result;
		}

		//------------------------------------------------compound assignment operators for Quaternion
		///quaternion + quaternion
		template<typename Type>
		inline TQuaternionExpression<TQuaternion<Type>, CQuaternionAddition, TQuaternion<Type> >
			operator + (const TQuaternion<Type>& lhs_, const TQuaternion<Type>& rhs_)
		{
			return TQuaternionExpression<TQuaternion<Type>, CQuaternionAddition, TQuaternion<Type> >(lhs_, rhs_);
		}

		///quaternion + expression
		template<typename Type, class Rhs>
		inline TQuaternionExpression<TQuaternion<Type>, CQuaternionAddition, Rhs>
			operator + (const TQuaternion<Type>& lhs_, const Rhs& rhs_)
		{
			return TQuaternionExpression<TQuaternion<Type>, CQuaternionAddition, Rhs>(lhs_, rhs_);
		}

		///expression + quaternion
		template<typename Type, class Lhs>
		inline TQuaternionExpression<Lhs, CQuaternionAddition, TQuaternion<Type> >
			operator + (const Lhs& lhs_, const TQuaternion<Type>& rhs_)
		{
			return TQuaternionExpression<Lhs, CQuaternionAddition, TQuaternion<Type> >(lhs_, rhs_);
		}

		///expression + expression
		template<class LLhs, typename LExp, class LRhs, class RLhs, typename RExp, class RRhs>
		inline TMatrixExpression<
			TQuaternionExpression<LLhs, LExp, LRhs>,
			CQuaternionAddition,
			TQuaternionExpression<RLhs, RExp, RRhs> >
		operator + (
			const TQuaternionExpression<LLhs, LExp, LRhs>& lhs_,
			const TQuaternionExpression<RLhs, RExp, RRhs>& rhs_)
		{
			return TMatrixExpression<
				TQuaternionExpression<LLhs, LExp, LRhs>,
				CQuaternionAddition,
				TQuaternionExpression<RLhs, RExp, RRhs> >(lhs_, rhs_);
		}

		///quaternion - quaternion
		template<typename Type>
		inline TQuaternionExpression<TQuaternion<Type>, CQuaternionSubtraction, TQuaternion<Type> >
			operator - (const TQuaternion<Type>& lhs_, const TQuaternion<Type>& rhs_)
		{
			return TQuaternionExpression<TQuaternion<Type>, CQuaternionSubtraction, TQuaternion<Type> >(lhs_, rhs_);
		}

		///quaternion - expression
		template<typename Type, class Rhs>
		inline TQuaternionExpression<TQuaternion<Type>, CQuaternionSubtraction, Rhs>
			operator - (const TQuaternion<Type>& lhs_, const Rhs& rhs_)
		{
			return TQuaternionExpression<TQuaternion<Type>, CQuaternionSubtraction, Rhs>(lhs_, rhs_);
		}

		///expression - quaternion
		template<typename Type, class Lhs>
		inline TQuaternionExpression<Lhs, CQuaternionSubtraction, TQuaternion<Type> >
			operator - (const Lhs& lhs_, const TQuaternion<Type>& rhs_)
		{
			return TQuaternionExpression<Lhs, CQuaternionSubtraction, TQuaternion<Type> >(lhs_, rhs_);
		}

		///expression - expression
		template<class LLhs, typename LExp, class LRhs, class RLhs, typename RExp, class RRhs>
		inline TMatrixExpression<
			TQuaternionExpression<LLhs, LExp, LRhs>,
			CQuaternionSubtraction,
			TQuaternionExpression<RLhs, RExp, RRhs> >
		operator - (
			const TQuaternionExpression<LLhs, LExp, LRhs>& lhs_,
			const TQuaternionExpression<RLhs, RExp, RRhs>& rhs_)
		{
			return TMatrixExpression<
				TQuaternionExpression<LLhs, LExp, LRhs>,
				CQuaternionSubtraction,
				TQuaternionExpression<RLhs, RExp, RRhs> >(lhs_, rhs_);
		}

		///quaternion * quaternion
		///
		///[w1, v1][w2, v2] = [w1w2-v1・v2, w1v2 + w2v1 + v2×v1]
		///
		template<typename Type>
		inline TQuaternion<Type>
			operator * (
				const TQuaternion<Type>& lhs_,
				const TQuaternion<Type>& rhs_) 
		{
			//create temporary object for reduce quaternion multiplication overhead
			TQuaternion<Type> tmp;

			tmp.m_x = lhs_.m_w*rhs_.m_x + lhs_.m_x*rhs_.m_w + lhs_.m_y*rhs_.m_z - lhs_.m_z*rhs_.m_y;
			tmp.m_y = lhs_.m_w*rhs_.m_y + lhs_.m_y*rhs_.m_w - lhs_.m_x*rhs_.m_z + lhs_.m_z*rhs_.m_x;
			tmp.m_z = lhs_.m_w*rhs_.m_z + lhs_.m_z*rhs_.m_w + lhs_.m_x*rhs_.m_y - lhs_.m_y*rhs_.m_x;
			tmp.m_w = lhs_.m_w*rhs_.m_w - lhs_.m_x*rhs_.m_x - lhs_.m_y*rhs_.m_y - lhs_.m_z*rhs_.m_z;
			return tmp;
		}

		///quaternion * expression
		template<typename Type, typename Lhs>
		inline TQuaternion<Type>
			operator * (
				const Lhs& lhs_,
				const TQuaternion<Type>& rhs_) 
		{
			//evaluate expression
			TQuaternion<Type> lhs(lhs_);

			//create temporary object for reduce quaternion multiplication overhead
			TQuaternion<Type> tmp;

			tmp.m_x = lhs.m_w*rhs_.m_x + lhs.m_x*rhs_.m_w + lhs.m_y*rhs_.m_z - lhs.m_z*rhs_.m_y;
			tmp.m_y = lhs.m_w*rhs_.m_y + lhs.m_y*rhs_.m_w - lhs.m_x*rhs_.m_z + lhs.m_z*rhs_.m_x;
			tmp.m_z = lhs.m_w*rhs_.m_z + lhs.m_z*rhs_.m_w + lhs.m_x*rhs_.m_y - lhs.m_y*rhs_.m_x;
			tmp.m_w = lhs.m_w*rhs_.m_w - lhs.m_x*rhs_.m_x - lhs.m_y*rhs_.m_y - lhs.m_z*rhs_.m_z;
			return tmp;
		}

		///expression * quaternion
		template<typename Type, typename Rhs>
		inline TQuaternion<Type>
			operator * (
				const TQuaternion<Type>& lhs_,
				const Rhs& rhs_) 
		{
			//evaluate expression
			TQuaternion<Type> rhs(rhs_);

			//create temporary object for reduce quaternion multiplication overhead
			TQuaternion<Type> tmp;

			tmp.m_x = lhs_.m_w*rhs.m_x + lhs_.m_x*rhs.m_w + lhs_.m_y*rhs.m_z - lhs_.m_z*rhs.m_y;
			tmp.m_y = lhs_.m_w*rhs.m_y + lhs_.m_y*rhs.m_w - lhs_.m_x*rhs.m_z + lhs_.m_z*rhs.m_x;
			tmp.m_z = lhs_.m_w*rhs.m_z + lhs_.m_z*rhs.m_w + lhs_.m_x*rhs.m_y - lhs_.m_y*rhs.m_x;
			tmp.m_w = lhs_.m_w*rhs.m_w - lhs_.m_x*rhs.m_x - lhs_.m_y*rhs.m_y - lhs_.m_z*rhs.m_z;
			return tmp;
		}

		///expression * expression
		template<class LLhs, typename LExp, class LRhs, class RLhs, typename RExp, class RRhs>
		inline TQuaternion<typename LLhs::value_type>
			operator * (
				const TQuaternionExpression<LLhs, LExp, LRhs>& lhs_,
				const TQuaternionExpression<RLhs, RExp, RRhs>& rhs_)
		{
			//create temporary object for reduce matrix multiplication overhead
			TQuaternion<LLhs::value_type> lhs(lhs_);
			TQuaternion<LLhs::value_type> rhs(rhs_);

			//create temporary object for reduce quaternion multiplication overhead
			TQuaternion<LLhs::value_type> tmp;

			tmp.m_x = lhs_.m_w*rhs_.m_x + lhs_.m_x*rhs_.m_w + lhs_.m_y*rhs_.m_z - lhs_.m_z*rhs_.m_y;
			tmp.m_y = lhs_.m_w*rhs_.m_y + lhs_.m_y*rhs_.m_w - lhs_.m_x*rhs_.m_z + lhs_.m_z*rhs_.m_x;
			tmp.m_z = lhs_.m_w*rhs_.m_z + lhs_.m_z*rhs_.m_w + lhs_.m_x*rhs_.m_y - lhs_.m_y*rhs_.m_x;
			tmp.m_w = lhs_.m_w*rhs_.m_w - lhs_.m_x*rhs_.m_x - lhs_.m_y*rhs_.m_y - lhs_.m_z*rhs_.m_z;
			return obj;
		}

		///quaternion * scalar
		template<typename T>
		inline TQuaternion<T> operator * (
			const TQuaternion<T>& lhs_,
			const typename TQuaternion<T>::value_type& rhs_)
		{	
			TQuaternion<T> obj(lhs_);
			obj.m_x *= rhs_;
			obj.m_y *= rhs_;
			obj.m_z *= rhs_;
			obj.m_w *= rhs_;
			return obj;
		}

		///scalar * quaternion
		template<typename T>
		inline TQuaternion<T> operator * (
			const typename TQuaternion<T>::value_type& lhs_,
			const TQuaternion<T>& rhs_)
		{	
			TQuaternion<T> obj(rhs_);
			obj.m_x *= lhs_;
			obj.m_y *= lhs_;
			obj.m_z *= lhs_;
			obj.m_w *= lhs_;
			return obj;
		}

		//------------------------------------------------combinations
		///vector * matrix
		///ただし、行列はDim×Dimの正方行列に限る
		///ETを使わず、一時オブジェクトに乗算の結果を保存しておいたおいたほうが
		///計算量が減ることがある。ここは最高速と平均速度のトレードオフになると思われる。
		///乗算は*=しか定義しない手もありえるか
		template<typename Type, int Row, int Colmun>
		inline TVector<Type, Colmun> 
			operator * (const TVector<Type, Row>& lhs_, const TMatrix<Type, Row, Colmun>& rhs_) 
		{
			TVector<Type, Colmun> obj;
			for (int i = 0; i < Colmun; ++i)
			{
				obj(i) = lhs_(0) * rhs_(0, i);
				for (int r = 1; r < Row; ++r)
					obj(i) += lhs_(r) * rhs_(r, i);
			}
			return obj;
		}

		///vector * matrix expression

		///vector expression * matrix

		///vector expression * matrix expression


	}
}










#endif