# MHD_solver

HLLD近似リーマン解法による磁気流体計算コードです。
現在以下をサンプルとして収録しています。

- Orszag-Tang渦問題 (Orszag and Tang, 1979)
- Kelvin-Helmholtz不安定 (Matsumoto and Hoshino, 2004)

計算方法等の詳細は[弊記事](https://qiita.com/ur_kinsk/items/1893602e2ee73060b207)をご覧ください。


## C++バージョンについて
コンパイルにはC++17以上が必要です。
また、`-fopenmp` オプションを有効にすると並列計算を行います。
Makefile等は未作成なので、以下のようにコンパイルしてください。

```
$ g++ -std=c++17 -O3 -fopenmp mhd.cc main.cc
```

実行すると `data` ディレクトリ内にアウトプットファイルが書き込まれるので、それを `plot_img.py` および `plot_mov.py` によって可視化できます。前者は `img` ディレクトリに画像を、後者は `mov` ディレクトリに動画を出力します。

## Fortranバージョンについて
使用方法はC++バージョンと同様です。コンパイルは以下のようにしてください。

```
$ gfortran -O3 -fopenmp mhd.f90 main.f90
```

## Pythonバージョンについて
計算終了後、その場で圧力をプロットします。
ファイル出力に関するコメントアウトを外せば、C++バージョンと同様 `data` ディレクトリにファイルが出力されるようになります。

## Juliaバージョンについて
使用方法はPythonバージョンと同様です。
パッケージとして `PyCall`, `CSV`, `DataFrames`, `Formatting` を事前にインストールしておいてください。
プロットには `matplotlib` を呼び出しています。
