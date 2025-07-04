import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import xml.etree.ElementTree as ET

class RaceLineEditApp:
    def __init__(self, root):
        self.root = root
        self.root.title("対話型走行ラインエディタ")
        self.root.geometry("1200x800")

        # 状態変数の初期化
        self.df = None; self.map_nodes = None; self.selected_item = None
        self.raceline_scatter = None; self.raceline_line = None # 走行ラインの「線」を管理
        self.current_annotation = None
        self.canvas = None; self.toolbar = None
        self.drag_data = {"index": None, "press_xy": None, "is_dragging": False}

        self.setup_left_panel()
        self.update_plot()
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def on_closing(self):
        plt.close('all')
        self.root.destroy()
    def setup_left_panel(self):
        """左側の操作パネルを構築する"""
        left_frame = ttk.Frame(self.root, padding="10"); left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        file_frame = ttk.LabelFrame(left_frame, text="ファイル操作", padding="10"); file_frame.pack(fill=tk.X, pady=5)
        ttk.Button(file_frame, text="走行ラインCSVを読み込む", command=self.load_csv).pack(fill=tk.X, pady=2)
        ttk.Button(file_frame, text="OSM地図を読み込む", command=self.load_osm_map).pack(fill=tk.X, pady=2)
        ttk.Button(file_frame, text="変更をCSVに保存", command=self.save_csv).pack(fill=tk.X, pady=2, padx=5)
        self.tree = ttk.Treeview(left_frame, columns=('Index', 'X', 'Y', 'Speed'), show='headings'); self.tree.heading('Index', text='点 #'); self.tree.heading('X', text='X座標'); self.tree.heading('Y', text='Y座標'); self.tree.heading('Speed', text='速度'); self.tree.column('Index', width=50); self.tree.column('X', width=100); self.tree.column('Y', width=100); self.tree.column('Speed', width=80); self.tree.pack(fill=tk.BOTH, expand=True, pady=10); self.tree.bind('<<TreeviewSelect>>', self.on_item_select)
        edit_frame = ttk.LabelFrame(left_frame, text="選択した点を編集", padding="10"); edit_frame.pack(fill=tk.X, pady=10)
        self.x_var, self.y_var, self.speed_var = tk.StringVar(), tk.StringVar(), tk.StringVar()
        ttk.Label(edit_frame, text="X座標:").grid(row=0, column=0, sticky=tk.W, pady=2); ttk.Entry(edit_frame, textvariable=self.x_var).grid(row=0, column=1, sticky=tk.EW)
        ttk.Label(edit_frame, text="Y座標:").grid(row=1, column=0, sticky=tk.W, pady=2); ttk.Entry(edit_frame, textvariable=self.y_var).grid(row=1, column=1, sticky=tk.EW)
        ttk.Label(edit_frame, text="速度:").grid(row=2, column=0, sticky=tk.W, pady=2); ttk.Entry(edit_frame, textvariable=self.speed_var).grid(row=2, column=1, sticky=tk.EW)
        ttk.Button(edit_frame, text="値を更新してプロット再描画", command=self.update_values_from_entry).grid(row=3, column=0, columnspan=2, pady=10)
        self.right_frame = ttk.Frame(self.root, padding="10"); self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    def update_plot(self):
        for widget in self.right_frame.winfo_children(): widget.destroy()
        fig, ax = plt.subplots(figsize=(8, 8)); self.current_annotation = None
        if self.map_nodes:
            map_x, map_y = zip(*self.map_nodes)
            ax.scatter(map_x, map_y, c='black', s=2, alpha=0.5, label='map')
        if self.df is not None and not self.df.empty:
            data = self.df[['x', 'y', 'speed']].to_numpy()
            self.raceline_scatter = ax.scatter(data[:, 0], data[:, 1], c=data[:, 2], cmap='jet', s=25, label='racing line(point)', picker=5, zorder=10)
            
            # ### イヤーッ！ これが【ウロボロス・ライン】のジツだ！ ###
            # プロット用のデータを作成し、点0の座標を末尾に追加する
            line_plot_data = np.vstack([data[:, :2], data[0, :2]])
            self.raceline_line, = ax.plot(line_plot_data[:, 0], line_plot_data[:, 1], color='cyan', linewidth=1.2, label='racing line (line)', zorder=5)

            cbar = fig.colorbar(self.raceline_scatter, ax=ax, shrink=0.8)
            cbar.set_label('V (m/s)', rotation=270, labelpad=20)
        ax.set_title(""); ax.set_xlabel("X (m)"); ax.set_ylabel("Y (m)")
        ax.grid(True); ax.set_aspect('equal', adjustable='box'); ax.legend()
        self.canvas = FigureCanvasTkAgg(fig, master=self.right_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.right_frame); self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
  
    def on_scroll(self, event):
        if event.inaxes is None or not hasattr(self, 'toolbar') or self.toolbar.mode: return
        ax = event.inaxes; cur_xlim = ax.get_xlim(); cur_ylim = ax.get_ylim(); xdata, ydata = event.xdata, event.ydata
        base_scale = 1.2
        if event.button == 'up': scale_factor = 1 / base_scale
        elif event.button == 'down': scale_factor = base_scale
        else: return
        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor; new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor
        rel_x = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0]); rel_y = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])
        ax.set_xlim([xdata - new_width * (1 - rel_x), xdata + new_width * rel_x]); ax.set_ylim([ydata - new_height * (1 - rel_y), ydata + new_height * rel_y])
        self.canvas.draw_idle()

    def on_press(self, event):
        if event.inaxes is None or self.raceline_scatter is None: return
        if self.toolbar.mode: return
        contains, attrd = self.raceline_scatter.contains(event)
        if contains:
            self.drag_data["index"] = attrd["ind"][0]
            self.drag_data["press_xy"] = (event.x, event.y)

    def on_motion(self, event):
        if self.drag_data["index"] is None or event.inaxes is None: return
        if not self.drag_data["is_dragging"]:
            px, py = self.drag_data["press_xy"]
            if (event.x - px)**2 + (event.y - py)**2 > 25:
                self.drag_data["is_dragging"] = True
                if self.current_annotation: self.current_annotation.remove(); self.current_annotation = None
        if self.drag_data["is_dragging"]:
            idx = self.drag_data["index"]
            offsets = self.raceline_scatter.get_offsets()
            offsets[idx] = [event.xdata, event.ydata]
            self.raceline_scatter.set_offsets(offsets)
            line_data = self.raceline_line.get_xydata()
            line_data[idx] = [event.xdata, event.ydata]
            # 始点または終点が動いた場合、閉じる線も更新
            if idx == 0: line_data[-1] = [event.xdata, event.ydata]
            if idx == len(self.df) -1: line_data[0] = [event.xdata, event.ydata]
            self.raceline_line.set_data(line_data[:, 0], line_data[:, 1])
            self.canvas.draw_idle()

    def on_release(self, event):
        if self.drag_data["index"] is None: return
        idx = self.drag_data["index"]
        if self.drag_data["is_dragging"]:
            self.df.loc[idx, 'x'] = event.xdata
            self.df.loc[idx, 'y'] = event.ydata
            self.populate_table(); self.highlight_and_focus_table_row(idx)
        else: # クリック操作
            if self.current_annotation: self.current_annotation.remove()
            ax = self.raceline_scatter.axes
            x_coord, y_coord = self.df.iloc[idx]['x'], self.df.iloc[idx]['y']
            self.current_annotation = ax.annotate(f"#{idx}", xy=(x_coord, y_coord), xytext=(15, 15), textcoords='offset points', bbox=dict(boxstyle="round,pad=0.4", fc="yellow", ec="black", lw=1, alpha=0.9), arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.1"))
            self.canvas.draw(); self.highlight_and_focus_table_row(idx)
        self.drag_data = {"index": None, "press_xy": None, "is_dragging": False}
    
    def highlight_and_focus_table_row(self, index_to_find):
        for iid in self.tree.get_children():
            if int(self.tree.item(iid, 'values')[0]) == index_to_find:
                self.tree.selection_set(iid); self.tree.focus(iid); self.tree.see(iid)
                self.on_item_select(None); break

    def update_values_from_entry(self):
        if self.selected_item is None or self.df is None: messagebox.showwarning("注意", "先にテーブルから更新したい点を選択してください。"); return
        try:
            idx = int(self.tree.item(self.selected_item, 'values')[0])
            new_x, new_y, new_speed = float(self.x_var.get()), float(self.y_var.get()), float(self.speed_var.get())
            self.df.loc[idx, 'x'], self.df.loc[idx, 'y'], self.df.loc[idx, 'speed'] = new_x, new_y, new_speed
            self.populate_table(); self.update_plot()
        except ValueError: messagebox.showerror("エラー", "X, Y, 速度には数値を入力してください。")
        except Exception as e: messagebox.showerror("エラー", f"予期せぬエラーが発生しました: {e}")

    def load_osm_map(self):
        file_path = filedialog.askopenfilename(filetypes=[("OSM/XML files", "*.osm *.txt *.xml")]);
        if not file_path: return
        try:
            tree = ET.parse(file_path); root = tree.getroot(); nodes = []
            for node in root.findall('node'):
                local_x, local_y = None, None
                for tag in node.findall('tag'):
                    if tag.get('k') == 'local_x': local_x = float(tag.get('v'))
                    if tag.get('k') == 'local_y': local_y = float(tag.get('v'))
                if local_x is not None and local_y is not None: nodes.append((local_x, local_y))
            self.map_nodes = nodes; messagebox.showinfo("成功", f"{len(self.map_nodes)}個の地図ノードを読み込みました。"); self.update_plot()
        except Exception as e: messagebox.showerror("エラー", f"地図ファイルの解析中にエラーが発生しました: {e}")

    # ### イヤーッ！ これが【円環の理】のジツだ！ ###
    def save_csv(self):
        """保存時に点0のコピーを末尾に加える"""
        if self.df is None: messagebox.showwarning("注意", "保存するデータがありません。"); return
        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")], title="名前を付けて保存", initialfile="raceline_modified.csv")
            if file_path:
                # 点0のコピーを末尾に加えた、保存専用のデータフレームを作成する
                df_to_save = pd.concat([self.df, self.df.iloc[[0]]], ignore_index=True)
                df_to_save.to_csv(file_path, index=False)
                messagebox.showinfo("成功", f"ファイルは正常に保存されました:\n{file_path}")
        except Exception as e: messagebox.showerror("エラー", f"ファイルの保存中にエラーが発生しました: {e}")

    def load_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")]);
        if not file_path: return
        self.df = pd.read_csv(file_path)
        # 始点と終点が同じ座標なら、終点を削除して円環を閉じる前の状態に戻す
        if len(self.df) > 1:
            last_idx = len(self.df) - 1
            if all(self.df.iloc[0] == self.df.iloc[last_idx]):
                self.df = self.df.iloc[:last_idx]
        rename_map = {'# x_m': 'x', 'x_m': 'x', 'y_m': 'y', 'vx_mps': 'speed'}; self.df.rename(columns=rename_map, inplace=True)
        if 'x' not in self.df.columns or 'y' not in self.df.columns or 'speed' not in self.df.columns: messagebox.showerror("エラー", "CSVに 'x', 'y', 'speed' のいずれかの列が含まれていません。"); self.df = None; return
        self.populate_table(); self.update_plot()

    def populate_table(self):
        if self.tree.get_children(): self.tree.delete(*self.tree.get_children())
        for index, row in self.df.iterrows():
            self.tree.insert("", "end", values=(index, f"{row['x']:.2f}", f"{row['y']:.2f}", f"{row['speed']:.2f}"))
    
    def on_item_select(self, event):
        if self.tree.selection():
            self.selected_item = self.tree.selection()[0]
            values = self.tree.item(self.selected_item, 'values')
            index, x, y, speed = values
            self.x_var.set(x); self.y_var.set(y); self.speed_var.set(speed)

if __name__ == "__main__":
    root = tk.Tk()
    app = RaceLineEditApp(root)
    root.mainloop()