rfmath
======

C++ math tool for my game project.


3Dグラフィックスの計算に必要な、ベクトル、クオータニオンの計算を行うライブラリです。
自作ゲーム向けなので、そこそこ正しく動けばいいやというスタンスで作成されています。
特殊なチューニング当はおこなっていませんが、逆に中身の処理は単調で理解しやすいかもしれません。だといいな。


■使い方

ライブラリを利用するにはmath/index.hをインクルードしてください。

■その他情報

index.h以外の個別ファイルの機能は以下


// 行列、ベクトル、クオータニオン、オイラー角

// 値を保持することに徹しており、特に難しい計算を直接扱うことはない。

vectortemplate.h

matrixtemplate.h

euleranglestemplate.h

quaterniontemplate.h


// sqrtやsin,conなど
math.h

// 行列の四則演算 式テンプレート

expressionoperators.h

expressiontemplate.h


// 内積とかのような四則演算のその上

linearalgebras.h


