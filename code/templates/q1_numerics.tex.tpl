\subsection{模型求解}

对径向区域作均匀离散，温度采用中心差分，水分采用\textbf{守恒通量离散}，使相邻控制体共用同一界面通量。圆心按轴对称极限处理，表面施加对流条件。空间离散后，采用\textbf{BDF 自适应积分}求解。

径向网格、控制体及热量和水分的界面通量关系见图~\ref{fig:radial-fvm-grid}。

\begin{figure}[H]
\centering
\includegraphics[width=0.94\textwidth]{figures/radial_fvm_grid.png}
\caption{圆柱药材的径向网格、控制体与界面通量示意图}
\label{fig:radial-fvm-grid}
\end{figure}

为分辨初期表面附近的含水率变化，逐级加密网格，最终采用 $@@N@@$ 个径向区间；相对和绝对时间容差分别为 $10^{-11}$、$10^{-13}$，最大时间步为 $5\,\mathrm{s}$。按题意每隔 $1\,\mathrm{s}$、每隔 $0.1\,\mathrm{cm}$ 输出，得到两组 $1801\times21$ 的状态值。\textbf{每秒输出不代表以一秒为积分步长}，内部步长由误差控制自适应确定。

\noindent\textbf{数值检验。}
保持时间容差一致，将径向区间数由 $@@PREV_N@@$ 增至 $@@N@@$，全部规定输出点的温度与含水率最大差异分别为 $@@SPACE_T@@\,{}^\circ\mathrm{C}$ 和 $@@SPACE_C@@\,\mathrm{kg/kg}$，\textbf{两级网格的规定结果表在四位小数下保持一致}。固定细网格并将两种时间容差同时收紧十倍后，两场最大差异分别为 $@@TIME_T@@\,{}^\circ\mathrm{C}$ 和 $@@TIME_C@@\,\mathrm{kg/kg}$。平均含水率变化与累计表面通量的守恒残差不超过 $@@MASS@@\,\mathrm{kg/kg}$。

另以均匀平衡解、圆柱对流导热的 Bessel 级数解析解及显式与隐式结果对照检查离散实现。上述检验支持数值求解的可靠性；端部与潜热简化的物理误差尚未量化。
