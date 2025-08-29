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
        self.root.title("raceline-visualizer for aaic2025")
        self.root.geometry("1200x800")

        # 状態変数の初期化
        self.df = None; self.map_nodes = None; self.selected_item = None
        self.raceline_scatter = None; self.raceline_line = None
        self.highlight_scatter = None
        self.current_annotation = None
        self.canvas = None; self.toolbar = None
        self.drag_data = {
            "indices": [],          # ドラッグ対象の点のインデックス(複数)
            "press_xdata_ydata": None, # ドラッグ開始時のマウスのデータ座標
            "initial_offsets": None, # ドラッグ開始時の各点の初期座標
            "is_dragging": False
        }
        self.last_selected_index = None

        self.history = [] 
        self.root.bind('<Control-z>', self.undo_last_action) 

        self.setup_left_panel()
        # self.root.update_idletasks()
        # self.update_plot()
        self.root.bind('<Configure>', self.on_first_draw)
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def on_closing(self):
        plt.close('all'); self.root.destroy()
        
    def setup_left_panel(self):
        """左側の操作パネルを構築する"""
        left_frame = ttk.Frame(self.root, padding="10"); left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        file_frame = ttk.LabelFrame(left_frame, text="ファイル操作", padding="10"); file_frame.pack(fill=tk.X, pady=5)
        ttk.Button(file_frame, text="走行ラインCSVを読み込む", command=self.load_csv).pack(fill=tk.X, pady=2)
        ttk.Button(file_frame, text="OSM地図を読み込む", command=self.load_osm_map).pack(fill=tk.X, pady=2)
        ttk.Button(file_frame, text="変更をCSVに保存", command=self.save_csv).pack(fill=tk.X, pady=2, padx=5)
        self.tree = ttk.Treeview(left_frame, columns=('Index', 'X', 'Y', 'Speed'), show='headings', selectmode='extended'); self.tree.heading('Index', text='点 #'); self.tree.heading('X', text='X座標'); self.tree.heading('Y', text='Y座標'); self.tree.heading('Speed', text='速度'); self.tree.column('Index', width=50); self.tree.column('X', width=100); self.tree.column('Y', width=100); self.tree.column('Speed', width=80); self.tree.pack(fill=tk.BOTH, expand=True, pady=10); self.tree.bind('<<TreeviewSelect>>', self.on_item_select)
        self.tree.bind('<Delete>', self.delete_points)
        edit_frame = ttk.LabelFrame(left_frame, text="選択した点を編集", padding="10"); edit_frame.pack(fill=tk.X, pady=10)
        self.x_var, self.y_var, self.speed_var = tk.StringVar(), tk.StringVar(), tk.StringVar()
        

        ttk.Label(edit_frame, text="X座標:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.x_entry = ttk.Entry(edit_frame, textvariable=self.x_var)
        self.x_entry.grid(row=0, column=1, sticky=tk.EW)
        
        ttk.Label(edit_frame, text="Y座標:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.y_entry = ttk.Entry(edit_frame, textvariable=self.y_var)
        self.y_entry.grid(row=1, column=1, sticky=tk.EW)


        ttk.Label(edit_frame, text="速度:").grid(row=2, column=0, sticky=tk.W, pady=2); ttk.Entry(edit_frame, textvariable=self.speed_var).grid(row=2, column=1, sticky=tk.EW)
        ttk.Button(edit_frame, text="値を更新してプロット再描画", command=self.update_values_from_entry).grid(row=3, column=0, columnspan=2, pady=5)
        ttk.Button(edit_frame, text="選択点の次に点を追加", command=self.add_point).grid(row=4, column=0, columnspan=2, pady=5)
        ttk.Button(edit_frame, text="選択した点を削除", command=self.delete_points).grid(row=5, column=0, columnspan=2, pady=5)
        
        interpolation_frame = ttk.LabelFrame(left_frame, text="補完機能", padding="10")
        interpolation_frame.pack(fill=tk.X, pady=10)
        ttk.Label(interpolation_frame, text="補間点数:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.interpolation_points_var = tk.StringVar(value="5")
        ttk.Entry(interpolation_frame, textvariable=self.interpolation_points_var).grid(row=0, column=1, sticky=tk.EW)
        ttk.Button(interpolation_frame, text="選択範囲を補間", command=self.interpolate_points).grid(row=1, column=0, columnspan=2, pady=5)

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
            line_plot_data = np.vstack([data[:, :2], data[0, :2]])
            self.raceline_line, = ax.plot(line_plot_data[:, 0], line_plot_data[:, 1], color='cyan', linewidth=1.2, label='racing line (line)', zorder=5)
            ####
            self.highlight_scatter, = ax.plot([], [], 'o', markersize=10, markerfacecolor='none', markeredgecolor='yellow', markeredgewidth=2, zorder=20)
            ####
            cbar = fig.colorbar(self.raceline_scatter, ax=ax, shrink=0.8); cbar.set_label('V (m/s)', rotation=270, labelpad=20)
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

    def _save_state_for_undo(self):
        """現在の状態を歴史の巻物に記録する"""
        if self.df is not None:
            self.history.append(self.df.copy())
            if len(self.history) > 100:
                self.history.pop(0)

    def undo_last_action(self, event=None):
        """時間を巻き戻し、一つ前の状態に戻す"""
        if self.history:
            self.df = self.history.pop()
            self.populate_table()
            self._redraw_plot_data_only()
            # Note: 巻き戻し後に選択状態を復元するのは複雑なので、ここでは行わない

    def on_first_draw(self, event=None):
        self.root.unbind('<Configure>')
        self.update_plot()

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
        if event.inaxes is None or self.raceline_scatter is None or self.toolbar.mode:
            return
        contains, attrd = self.raceline_scatter.contains(event)
        if contains:
            clicked_idx = attrd["ind"][0]

            # Shiftキーによる範囲選択
            if event.key == 'shift' and self.last_selected_index is not None:
                start_idx, end_idx = self.last_selected_index, clicked_idx
                num_points = len(self.df)
                dist_forward = (end_idx - start_idx + num_points) % num_points
                dist_backward = (start_idx - end_idx + num_points) % num_points
                indices = []
                if dist_forward <= dist_backward:
                    i = start_idx
                    while i != end_idx:
                        indices.append(i)
                        i = (i + 1) % num_points
                else:
                    i = start_idx
                    while i != end_idx:
                        indices.append(i)
                        i = (i - 1 + num_points) % num_points
                indices.append(end_idx)
                self.tree.selection_set([str(i) for i in indices])
                self.on_item_select(None) 
                return 
            
            self._save_state_for_undo()
            selected_iids = self.tree.selection()
            selected_indices = [int(iid) for iid in selected_iids]
            
            if clicked_idx in selected_indices:
                self.drag_data["indices"] = selected_indices
                self.drag_data["initial_offsets"] = self.df.loc[selected_indices, ['x', 'y']].to_numpy()
            else:
                self.drag_data["indices"] = [clicked_idx]
                self.drag_data["initial_offsets"] = self.df.loc[[clicked_idx], ['x', 'y']].to_numpy()
                self.highlight_and_focus_table_row(clicked_idx)

            self.drag_data["press_xdata_ydata"] = (event.xdata, event.ydata)
            self.last_selected_index = clicked_idx

    def on_motion(self, event):
        if not self.drag_data["indices"] or event.inaxes is None:
            return
        
        if not self.drag_data["is_dragging"]:
            if self.drag_data["press_xdata_ydata"]:
                dx = event.xdata - self.drag_data["press_xdata_ydata"][0]
                dy = event.ydata - self.drag_data["press_xdata_ydata"][1]
                if dx**2 + dy**2 > 1e-3: 
                    self.drag_data["is_dragging"] = True
                    if self.current_annotation:
                        self.current_annotation.remove()
                        self.current_annotation = None

        if self.drag_data["is_dragging"]:
            start_x, start_y = self.drag_data["press_xdata_ydata"]
            dx = event.xdata - start_x
            dy = event.ydata - start_y
            
            all_offsets = self.raceline_scatter.get_offsets()
            line_data = self.raceline_line.get_xydata()

            for i, idx in enumerate(self.drag_data["indices"]):
                initial_offset = self.drag_data["initial_offsets"][i]
                new_pos = [initial_offset[0] + dx, initial_offset[1] + dy]
                all_offsets[idx] = new_pos
                line_data[idx] = new_pos
                if idx == 0 and 0 in self.drag_data["indices"]:
                    line_data[-1] = new_pos
            
            self.raceline_scatter.set_offsets(all_offsets)
            self.raceline_line.set_data(line_data[:, 0], line_data[:, 1])

            highlight_indices = self.drag_data["indices"]
            highlight_coords = all_offsets[highlight_indices] 
            if highlight_coords.size > 0:
                self.highlight_scatter.set_data(highlight_coords[:, 0], highlight_coords[:, 1])
            else:
                self.highlight_scatter.set_data([], [])
            self.canvas.draw_idle()
    def on_release(self, event):
        if not self.drag_data["indices"]:
            return
            
        if self.drag_data["is_dragging"]:
            start_x, start_y = self.drag_data["press_xdata_ydata"]
            dx = event.xdata - start_x if event.xdata is not None else 0
            dy = event.ydata - start_y if event.ydata is not None else 0

            all_offsets = self.raceline_scatter.get_offsets()
            line_data = self.raceline_line.get_xydata()

            for i, idx in enumerate(self.drag_data["indices"]):
                initial_offset = self.drag_data["initial_offsets"][i]
                
                new_x = initial_offset[0] + dx
                new_y = initial_offset[1] + dy

                self.df.loc[idx, 'x'] = new_x
                self.df.loc[idx, 'y'] = new_y
                all_offsets[idx] = [new_x, new_y]
                line_data[idx] = [new_x, new_y]
                # 閉路の場合の始点/終点処理
                if idx == 0 and 0 in self.drag_data["indices"]:
                    line_data[-1] = [new_x, new_y]
            self.raceline_scatter.set_offsets(all_offsets)
            self.raceline_line.set_data(line_data[:, 0], line_data[:, 1])
            self.populate_table()
            self.tree.selection_set([str(i) for i in self.drag_data["indices"]])
            self.update_highlight() 

        else:
            if self.history:
                self.history.pop() 
            
            clicked_idx = self.drag_data["indices"][0]
            if self.current_annotation:
                self.current_annotation.remove()
            
            ax = self.raceline_scatter.axes
            point_data = self.df.iloc[clicked_idx]
            self.current_annotation = ax.annotate(f"#{clicked_idx}", 
                xy=(point_data['x'], point_data['y']), 
                xytext=(15, 15), textcoords='offset points',
                bbox=dict(boxstyle="round,pad=0.4", fc="yellow", ec="black", lw=1, alpha=0.9),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.1"))
            self.canvas.draw()
            self.highlight_and_focus_table_row(clicked_idx)
        
        # ドラッグ情報をリセット
        self.drag_data = {"indices": [], "press_xdata_ydata": None, "initial_offsets": None, "is_dragging": False}

    
    def highlight_and_focus_table_row(self, index_to_find):
        str_index = str(index_to_find)
        self.tree.selection_set(str_index)
        self.tree.focus(str_index)
        self.tree.see(str_index)
        self.on_item_select(None)

    def update_values_from_entry(self):
        selected_iids = self.tree.selection()
        if not selected_iids: messagebox.showwarning("注意", "先にテーブルから更新したい点を選択してください。"); return
        try:
            self._save_state_for_undo()
            new_speed = float(self.speed_var.get())
            for iid in selected_iids:
                idx = int(iid)
                self.df.loc[idx, 'speed'] = new_speed
            self.populate_table()
            self._redraw_plot_data_only() 
            self.tree.selection_set(selected_iids) 
        except ValueError: messagebox.showerror("エラー", "速度には数値を入力してください。")
        except Exception as e: messagebox.showerror("エラー", f"予期せぬエラーが発生しました: {e}")

    def add_point(self):
        selected_iids = self.tree.selection()
        if not selected_iids:
            messagebox.showwarning("注意", "点を追加したい場所の、直前の点を選択してください。")
            return

        self._save_state_for_undo()
        
        selected_indices = sorted([int(iid) for iid in selected_iids], reverse=True)

        for idx in selected_indices:
            p1 = self.df.iloc[idx]
            
            next_idx = (idx + 1) % len(self.df)
            p2 = self.df.iloc[next_idx]

            new_x = (p1['x'] + p2['x']) / 2
            new_y = (p1['y'] + p2['y']) / 2
            new_speed = (p1['speed'] + p2['speed']) / 2
            
            new_point_df = pd.DataFrame([{'x': new_x, 'y': new_y, 'speed': new_speed}])

            df_part1 = self.df.iloc[:idx + 1]
            df_part2 = self.df.iloc[idx + 1:]
            self.df = pd.concat([df_part1, new_point_df, df_part2]).reset_index(drop=True)

        self.populate_table()
        self._redraw_plot_data_only()

    def delete_points(self, event=None):
        selected_iids = self.tree.selection()
        if not selected_iids:
            messagebox.showwarning("注意", "削除したい点をテーブルで選択してください。")
            return

        if messagebox.askyesno("確認", f"{len(selected_iids)}個の点を選択しています。本当に削除しますか？"):
            self._save_state_for_undo()
            
            indices_to_delete = sorted([int(iid) for iid in selected_iids], reverse=True)
            
            self.df.drop(indices_to_delete, inplace=True)
            self.df.reset_index(drop=True, inplace=True)
            
            self.populate_table()
            self._redraw_plot_data_only()

    def catmull_rom_spline(self, P0, P1, P2, P3, num_points):
        """Catmull-Romスプライン補間点を計算する"""
        alpha = 0.5  # Catmull-Romスプラインの種類を定義（0.5でCentripetal）
        
        def tj(ti, Pi, Pj):
            xi, yi = Pi
            xj, yj = Pj
            return ((xj - xi)**2 + (yj - yi)**2)**(alpha/2) + ti

        t0 = 0
        t1 = tj(t0, P0, P1)
        t2 = tj(t1, P1, P2)
        t3 = tj(t2, P2, P3)

        t = np.linspace(t1, t2, num_points + 2)[1:-1]
        
        # Handle cases where points are identical to avoid division by zero
        t1_t0 = (t1 - t0) if (t1 - t0) != 0 else 1e-6
        t2_t1 = (t2 - t1) if (t2 - t1) != 0 else 1e-6
        t3_t2 = (t3 - t2) if (t3 - t2) != 0 else 1e-6
        t2_t0 = (t2 - t0) if (t2 - t0) != 0 else 1e-6
        t3_t1 = (t3 - t1) if (t3 - t1) != 0 else 1e-6

        A1 = ((t1 - t)[:, np.newaxis] / t1_t0) * P0 + ((t - t0)[:, np.newaxis] / t1_t0) * P1
        A2 = ((t2 - t)[:, np.newaxis] / t2_t1) * P1 + ((t - t1)[:, np.newaxis] / t2_t1) * P2
        A3 = ((t3 - t)[:, np.newaxis] / t3_t2) * P2 + ((t - t2)[:, np.newaxis] / t3_t2) * P3
        
        B1 = ((t2 - t)[:, np.newaxis] / t2_t0) * A1 + ((t - t0)[:, np.newaxis] / t2_t0) * A2
        B2 = ((t3 - t)[:, np.newaxis] / t3_t1) * A2 + ((t - t1)[:, np.newaxis] / t3_t1) * A3

        C = ((t2 - t)[:, np.newaxis] / t2_t1) * B1 + ((t - t1)[:, np.newaxis] / t2_t1) * B2
        return C

    def interpolate_points(self):
        selected_iids = self.tree.selection()
        if len(selected_iids) < 2:
            messagebox.showwarning("注意", "補間するには、連続した2つ以上の点を選択してください。")
            return
        
        try:
            num_inserted_points = int(self.interpolation_points_var.get())
            if num_inserted_points < 1:
                messagebox.showerror("エラー", "補間点数は1以上で指定してください。")
                return
        except ValueError:
            messagebox.showerror("エラー", "補間点数には整数を入力してください。")
            return

        selected_indices = sorted([int(iid) for iid in selected_iids])

        # 選択が連続しているかチェック
        for i in range(len(selected_indices) - 1):
            if selected_indices[i+1] - selected_indices[i] != 1:
                messagebox.showwarning("注意", "補間するには、連続した点を選択してください。")
                return
        
        self._save_state_for_undo()
        
        all_new_points = []
        
        # Iterate through pairs of selected points to create segments
        for i in range(len(selected_indices) - 1):
            p1_series = self.df.iloc[selected_indices[i]]
            p2_series = self.df.iloc[selected_indices[i+1]]

            # Add the starting point of the segment as a dictionary
            if i == 0:
                all_new_points.append(p1_series.to_dict())

            # Interpolate all numeric columns
            for j in range(1, num_inserted_points + 1):
                t = j / (num_inserted_points + 1)
                new_point = {}
                for col in self.df.columns:
                    if pd.api.types.is_numeric_dtype(self.df[col]):
                        new_point[col] = p1_series[col] + (p2_series[col] - p1_series[col]) * t
                    else:
                        new_point[col] = p1_series[col]  # For non-numeric columns, copy from the first point
                all_new_points.append(new_point)
            
            # Add the ending point of the segment as a dictionary
            all_new_points.append(p2_series.to_dict())

        # Create a DataFrame from the new points
        new_df = pd.DataFrame(all_new_points).drop_duplicates().reset_index(drop=True)

        # Replace the old segment with the new interpolated segment in the original DataFrame
        df_part1 = self.df.iloc[:selected_indices[0]]
        df_part2 = self.df.iloc[selected_indices[-1] + 1:]

        self.df = pd.concat([df_part1, new_df, df_part2]).reset_index(drop=True)

        self.populate_table()
        self._redraw_plot_data_only()

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

    def save_csv(self):
        if self.df is None: messagebox.showwarning("注意", "保存するデータがありません。"); return
        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")], title="名前を付けて保存", initialfile="raceline_modified.csv")
            if file_path:
                df_to_save = pd.concat([self.df, self.df.iloc[[0]]], ignore_index=True)
                df_to_save.to_csv(file_path, index=False)
                messagebox.showinfo("成功", f"ファイルは正常に保存されました:\n{file_path}")
        except Exception as e: messagebox.showerror("エラー", f"ファイルの保存中にエラーが発生しました: {e}")

    def load_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")]);
        if not file_path: return
        try:
            self.df = pd.read_csv(file_path)
            if len(self.df) > 1:
                last_idx = len(self.df) - 1
                if all(self.df.iloc[0].round(5) == self.df.iloc[last_idx].round(5)):
                    self.df = self.df.iloc[:last_idx]
            rename_map = {'# x_m': 'x', 'x_m': 'x', 'y_m': 'y', 'vx_mps': 'speed'}; self.df.rename(columns=rename_map, inplace=True)
            if 'x' not in self.df.columns or 'y' not in self.df.columns or 'speed' not in self.df.columns: messagebox.showerror("エラー", "CSVに 'x', 'y', 'speed' のいずれかの列が含まれていません。"); self.df = None; return
            self.history = []
            self.populate_table(); self.update_plot()
        except pd.errors.ParserError:
            messagebox.showerror("エラー", "CSVファイルの読み込みに失敗しました。フォーマットを確認してください。")
        except Exception as e:
            messagebox.showerror("エラー", f"予期せぬエラーが発生しました: {e}")

    def populate_table(self):
        if self.tree.get_children(): self.tree.delete(*self.tree.get_children())
        if self.df is None: return
        for index, row in self.df.iterrows():
            self.tree.insert("", "end", iid=str(index), values=(index, f"{row['x']:.2f}", f"{row['y']:.2f}", f"{row['speed']:.2f}"))
            
    def update_highlight(self):
        if self.highlight_scatter is None: return
        selected_indices = [int(i) for i in self.tree.selection()]
        
        if selected_indices and self.df is not None:
            highlight_data = self.df.iloc[selected_indices][['x', 'y']].to_numpy()
            self.highlight_scatter.set_data(highlight_data[:, 0], highlight_data[:, 1])
        else:
            self.highlight_scatter.set_data([], [])   
        self.canvas.draw_idle()

    def on_item_select(self, event):
        selected_iids = self.tree.selection()
        num_selected = len(selected_iids)

        if num_selected > 1:
            # 複数選択時：座標をロックし、表示を明確化
            self.x_entry.config(state='disabled')
            self.y_entry.config(state='disabled')
            self.x_var.set("---")
            self.y_var.set("---")
        else:
            # 単数または選択なし：ロックを解除
            self.x_entry.config(state='normal')
            self.y_entry.config(state='normal')

        if selected_iids:
            # 選択がある場合、フォーカスされている行の値を表示する
            last_selected_iid = self.tree.focus()
            if not last_selected_iid: last_selected_iid = selected_iids[0]
            self.selected_item = last_selected_iid
            try:
                values = self.tree.item(self.selected_item, 'values')
                if values:
                    _, x, y, speed = values
                    # 単数選択の時だけX,Yをセットする
                    if num_selected == 1:
                        self.x_var.set(x)
                        self.y_var.set(y)
                    # 速度は常にフォーカスされている点の値を表示
                    self.speed_var.set(speed)
            except tk.TclError:
                pass
            self.update_highlight()
        else:
            # 選択がなくなった場合、テキストボックスを空にする
            self.x_var.set("")
            self.y_var.set("")
            self.speed_var.set("")
            self.update_highlight()
    def _redraw_plot_data_only(self):
        if self.df is None or not hasattr(self, 'raceline_scatter') or self.raceline_scatter is None:
            self.update_plot()
            return
        data = self.df[['x', 'y', 'speed']].to_numpy()
        if data.shape[0] > 0:
            new_min_speed = data[:, 2].min()
            new_max_speed = data[:, 2].max()
            self.raceline_scatter.norm.vmin = new_min_speed
            self.raceline_scatter.norm.vmax = new_max_speed
        self.raceline_scatter.set_offsets(data[:, :2])
        self.raceline_scatter.set_array(data[:, 2]) 
        line_plot_data = np.vstack([data[:, :2], data[0, :2]])
        self.raceline_line.set_data(line_plot_data[:, 0], line_plot_data[:, 1])
        ax = self.raceline_scatter.axes
        ax.relim()
        ax.autoscale_view()
        self.update_highlight()
if __name__ == "__main__":
    root = tk.Tk()
    app = RaceLineEditApp(root)
    root.mainloop()
