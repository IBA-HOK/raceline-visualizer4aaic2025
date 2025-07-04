# raceline-visualizer for aaic2025 

このツールは、自動運転AIチャレンジ2025用にMMTTチームが開発した、CSV形式の走行ラインデータを視覚的に表示し、対話的に編集するためのGUIアプリケーションです。

## 機能概要

本ツールは、CSV形式の走行ラインを対話的に編集するためのGUIアプリケーションです。主な機能として、CSV/OSMファイルの読み込み、速度に応じたラインの可視化、プロット上の点のドラッグ＆ドロップ編集、テーブルからの数値編集、マウスやツールバーによるプロット操作、編集結果の保存機能などを搭載しています。

---

## 環境構築 (Installation)

**必要なライブラリ:**

* Python 3.x
* pandas
* NumPy
* Matplotlib

以下のコマンドで、必要なライブラリをインストールしてください。

```bash
pip install pandas numpy matplotlib
```

**Ubuntu / Debian系OSでの注意点:**
お使いの環境によっては、`matplotlib`が`tkinter`の画像機能を正しく利用するために、追加のパッケージが必要になる場合があります。`ImportError: cannot import name 'ImageTk' from 'PIL'` のようなエラーが出た場合は、以下のコマンドを実行してください。

```bash
sudo apt-get install python3-pil.imagetk
```

---

## 使い方 (Usage)

1.  **起動**:
    ターミナルで以下のコマンドを実行して、アプリケーションを起動します。
    ```bash
    python3 plot.py
    ```
    (`plot.py`は、スクリプトを保存したファイル名に置き換えてください)

2.  **ファイルの読み込み**:
    * **「走行ラインCSVを読み込む」**: 編集したい走行ラインのCSVファイルを開きます。
    * **「OSM地図を読み込む」**: 背景として表示したい地図ファイル（`.osm`または`.txt`形式）を開きます。

3.  **走行ラインの編集**:
    * **ドラッグで編集**: グラフ上の点をマウスで掴んで、好きな位置へ移動させます。
    * **数値で編集**: 左のリストから編集したい点を選択し、その下に表示される「X座標」「Y座標」「速度」のボックスに直接数値を入力して「値を更新してプロット再描画」ボタンを押します。

4.  **保存**:
    * **「変更をCSVに保存」**: 編集後の走行ラインデータを新しいCSVファイルとして保存します。このとき、自動的に始点（点0）のデータが最終行に追加され、閉じたコースとして保存されます。

---

## ライセンス (License)

このプロジェクトは、以下のいずれかのライセンスを選択して利用できるデュアルライセンスです。

This project is dual-licensed, allowing you to choose between one of the following licenses.

### MIT License

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

### Creative Commons Zero v1.0 Universal (CC0)

あるいは、このプロジェクトはCC0の下でも利用可能です。これにより、著作権法の下で可能な限り、全ての権利を放棄します。あなたは、私たちの許可を得ることなく、商業目的であっても、このソフトウェアをコピー、変更、配布、実行することができます。
