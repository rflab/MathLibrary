//2010/04/18
//rflab.

#ifndef _RF_MATH_H_
#define _RF_MATH_H_

#include <cmath>

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
		///定数定義
		namespace definition
		{
			static const float RF_PI       = 3.14159265f; ///<π
			static const float RF_2PI      = 2*RF_PI;     ///<２π
			static const float RF_PI_OVER2 = RF_PI/2.0f;  ///<π/２
			static const float RF_OVER_PI  = 1.0f/RF_PI;  ///<１/π
			static const float RF_OVER_2PI = 1.0f/RF_2PI; ///<１/２π
			static const float RF_EPSILON  = 1e-4f;       ///<微小値
			static const int RF_UINT_MAX   = 0xffffffff;    ///<符号無し32ビット整数max値
		};

		///trigonometry
		class CTrigonometry
		{
		private:
			typedef enum{
				RESOLUTION = 0x400, //0x1ff+1//360度以上の回転はマスクで巡回できる値にしない解けない
				MASK       = 0x03ff,
				HALF       = RESOLUTION>>2
			}eDEFINE;

			//2/RF_PI*512
			const float RADTOBIT;
		private:
			float	m_sin[RESOLUTION];	//cosのテーブルをあらかじめ持っている

		private:
			CTrigonometry()
			:RADTOBIT(RESOLUTION/2/definition::RF_PI)
			{
				//make sin look up table
				for(int i=0; i<RESOLUTION; i++)
				{
					//2π/512した値を16ビット分整数値に押し上げる
					m_sin[i] = (std::sin)(i*2*definition::RF_PI/RESOLUTION);
				}
			}
		public:
			static CTrigonometry& Instance()
			{
				static CTrigonometry obj;
				return obj;
			}


			//本当はこのクラスは派生しないはずなのでvirtualは要らないはず
			virtual ~CTrigonometry(){}


			float Sin(float rad_) const 
			{
				return m_sin[((int)(rad_ * RADTOBIT)) & MASK];
			}
			
			float Cos(float rad_) const 
			{
				return m_sin[((int)((rad_ * RADTOBIT)+HALF)) & MASK];
			}

			void SinCos(float rad_, float* pSin_, float* pCos_) const 
			{
				*pSin_ = m_sin[((int)(rad_ * RADTOBIT)) & MASK];
				*pCos_ = m_sin[((int)((rad_ * RADTOBIT)+HALF)) & MASK];
			}

			//未実装orz
			static float Tan(float rad_) //const 
			{
				return (std::tan)(rad_);
			}

			static float ASin(float sin_) //const 
			{
				return asin(sin_);
			}

			static float ACos(float cos_) //const 
			{
				return acos(cos_);
			}


			static float ATan2(float y_, float x_) //const 
			{
				return atan2(y_, x_);
			}
		};

		///logarithm
		class CLogarithm
		{
		private:
			CLogarithm(){}
			virtual ~CLogarithm(){}

		public:
			static float Log(float x_) //const 
			{
				return (std::log)(x_);
			}

			static float Log10(float x_) //const 
			{
				return (std::log10)(x_);
			}
		};

		//-------------------------------------------------------randam
		class CRandom
		{
		private:
			mutable unsigned int m_seed;
			mutable unsigned int y;
			mutable unsigned int z;
			mutable unsigned int w;

			CRandom()
				:m_seed(0),
				 y(362436069), //10101100110100101010111100101 
				 z(521288629), //11111000100100011101110110101 
				 w(88675123)   //00101010010010001001100110011 
			{
			}

		public:
			virtual ~CRandom(){}

			static CRandom& Instance()	//簡易シングルトン
			{
				static CRandom obj;
				return obj;
			}

			void SetSeed(unsigned int seed_)
			{
				m_seed = seed_;
			}

			//Xorshiftアルゴリズム。理解してないorz
			unsigned int Rand() const
			{
				#if 1
					//unsigned int m_seed = 0;
					unsigned long t;

					t=m_seed^(m_seed<<11);
					m_seed=y;
					y=z;
					z=w;
					w=(w^(w>>19))^(t^(t>>8));
					return w;
				#else
					return rand();
				#endif
			}

			float RandUnit() const
			{
				return (float)Rand()/(float)definition::RF_UINT_MAX;
			}

			float RandUnitSigned() const
			{
				//気にする必要はないと思うが、処理系が変わると
				//unsigned int → intのキャストでマイナスにならないことがある？
				return (float)(Rand()/(float)definition::RF_UINT_MAX-0.5f)*2.0f;
			}
		};
		
		///square root
		class CSquareRoot
		{
			CSquareRoot(){}
			~CSquareRoot(){}
		public:

			///平方根
			static float Sqrt(float x_)// const
			{  
				#if 0		
					///平方根を求めるには　x2=c　のxを数値的に求めることができれば良いので
					///ニュートン法を使って求めることが出来ます。
					///漸化式は　xn = xn - (xn2 - c)/2xn

					//ニュートン法初期点の決定
					//see http://www15.plala.or.jp/hrak/contents/programdiary.html
					unsigned int i = *(unsigned int*) &x_;
					// adjust bias
					i  += 127 << 23;
					// approximation of square root
					i >>= 1;
					t = *(float*) &i;

					//ニュートン法初期点の決定
					if(x_<=0.0f)
						return 0.0f;


					//ニュートン法の実施
					float prev;
					do 
					{
						prev = t;
						t = (x_ / t + t) * 0.5f;
					}while (t < prev);

					return prev;
				#elif 1
					//高速化版
					//平方根の逆数に平方を掛ける
					return ReciprocalSqrt(x_)*x_;
				#else 
					//標準関数版
					return (std::sqrt)(x_);
				#endif
			}

			///1.0/sqrt(x)の高速化
			///f(X) = (1/x)^2 - c
			static float ReciprocalSqrt(float x_)// const 
			{
				#if 0
					//なぞのトリック
					//場合によって遅い
					unsigned int iDat;
					double dx,dy;

					iDat = ((0x7f000000 + 0x3f800000) - (*((unsigned int *)&x_))) >> 1;
					dy = *((float *)&iDat);

					dx = x_;
					dy = dy * 1.47 - dx * 0.47 * dy * dy * dy;
					dy = (3.0 - dy * dy * dx) * 0.5 * dy;

					return (float)dy;
				#elif 0
					//windowsOSへの最適化らしい
					unsigned int iDat;
					double dx,dy;

					if(x_<(1e-5)) return 0.0f;

					if(x_>(1e+5))
					{

						//少数部分をカット
						#if _WINDOWS
							x_ = (float)floor(x_);
						#else
							x_ = (float)((int)x_);
						#endif

						x_ = (float)((double)x_ * 0.00390625);  // x_=x_/256.0

						iDat = ((0x7f000000 + 0x3f800000) - (*((unsigned int *)&x_))) >> 1;
						dy = *((float *)&iDat);

						dx = x_;
						dy = dy * 1.47 - dx * 0.47 * dy * dy * dy;
						dy = (3.0 - dy * dy * dx) * 0.03125 * dy;
					}
					else
					{
						iDat = ((0x7f000000 + 0x3f800000) - (*((unsigned int *)&x_))) >> 1;
						dy = *((float *)&iDat);

						dx = x_;
						dy = dy * 1.47 - dx * 0.47 * dy * dy * dy;
						dy = (3.0 - dy * dy * dx) * 0.5 * dy;
					}

					return (float)dy;
				#elif 1
					//Quake3の高速逆平方根トリック
					float halfX = 0.5f*x_;
					int i = *(int*)&x_;
					i = 0x5f3759df - (i>>1);
					x_ = *(float*)&i;
					return x_*(1.5f - halfX*x_*x_);
				#else
					//標準関数版
					return 1.0f/sqrt(x_);
				#endif
			}
		};

		///グローバル関数
		namespace function
		{
			//-----------------------------trigonometry
			inline float Sin(float rad_){return CTrigonometry::Instance().Sin(rad_);}
			inline float Cos(float rad_){return CTrigonometry::Instance().Cos(rad_);}
			inline void SinCos(float rad_, float* pSin_, float* pCos_)
				{return CTrigonometry::Instance().SinCos(rad_, pSin_, pCos_);}
			inline float Tan(float rad_){return CTrigonometry::Tan(rad_);}
			inline float ASin(float sin_){return CTrigonometry::ASin(sin_);}
			inline float ACos(float cos_){return CTrigonometry::ACos(cos_);}
			inline float ATan2(float y_, float x_){return CTrigonometry::ATan2(y_, x_);}
	
			//-----------------------------logarithm
			inline float Log(float x_){return CLogarithm::Log(x_);}
			inline float Log10(float x_){return CLogarithm::Log10(x_);}

			//-----------------------------randam
			inline unsigned int RandUnit(){return CRandom::Instance().Rand();}
			inline float RandUnitSigned(){return CRandom::Instance().RandUnit();}

			//-----------------------------square root
			inline float Sqrt(float x_){return CSquareRoot::Sqrt(x_);}
			inline float RSqrt(float x_){return CSquareRoot::ReciprocalSqrt(x_);}

			//-----------------------------others

			///最小、最大の上限にクランプ
			template<typename T>
			inline void Clamp(T* target_, const T& min_, const T& max_)	
			{
				if (*target_ < min_)
					*target_ = min_;
				else if (*target_ > max_)
					*target_ = max_;
			}

			///絶対値関数
			template<typename T>
			inline T Abs(T value_)
			{
				return std::abs(value_);
			}
			template<>
			inline float Abs(float value_)
			{
				return std::fabs(value_);
			}

			///床関数
			template<typename T>
			inline T Floor(T value_)
			{
				return std::floor(value_);
			}

			///天井関数
			template<typename T>
			inline T Ceil(T value_)
			{
				return std::ceil(value_);
			}

			///min関数
			template<typename T>
			inline T Min(T _a, T _b)
			{
				return (std::min)(_a, _b);
			}

			///max関数
			template<typename T>
			inline T Max(T _a, T _b)
			{
				return (std::max)(_a, _b);
			}

			///最小最大の範囲で周回させた値を得る
			template<typename T>
			inline void Cycle(T* target_, const T& min_, const T& max_)	
			{
				//最小値にオフセットし
				//レンジで割った小数部分にレンジをかけてレンジ内にラップする
				//誤差が大きそうだけど3Dmathのソースコードでこうやってた
				T range = max_-min_;
				*target_ -= min_;
				*target_ -= Floor(target_/range)*range; 
				*target_ -= min_;
			}
		}
	}
}


#endif