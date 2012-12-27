//2010/04/18
//rflab.

#ifndef _RF_SCOPE_TIME_H_
#define _RF_SCOPE_TIME_H_

#include <time.h>
#include <string.h>
#include <iostream>
#include <sstream>

#include <stdio.h>

//debug
#undef LOCALPRINTF
//#define LOCALPRINTF(...) printf(__VA_ARGS__)
#define LOCALPRINTF(...)

///rock field!
namespace rf
{
	///時間計測クラス                       <br>
	///インスタンス生成からの時間を計測する <br>
	///基本的にデバッグ用                   <br>
	class CClock
	{
		clock_t m_start;

	public:
		///LOCALPRINTFを有効にした場合、デストラクタで時間を表示
		~CClock()
		{
			LOCALPRINTF("at destructor %f[sec.]\n", (float)(clock()-m_start)/(float)CLOCKS_PER_SEC);
		}

		CClock()
		{
			Reset();
		}

		///現時点の時間を文字列に変換
		operator std::string()
		{
			clock_t end = clock();
			std::stringstream ss;
			ss << (float)(end-m_start)/(float)CLOCKS_PER_SEC << "[sec.]";
			return ss.str();
		}

		///現時点の時間をfloat値に変換
		operator float ()
		{
			return (float)(clock()-m_start)/(float)CLOCKS_PER_SEC;
		}
		
		///再計測開始
		void Reset()
		{
			m_start = clock();
		}
	};
}

#endif
