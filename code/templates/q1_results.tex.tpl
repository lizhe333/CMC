\Needspace{6\baselineskip}
\subsection{结果分析}

\noindent\textbf{（1）规定时空点结果}

规定时刻和径向位置的温度、干基含水率分别见表~\ref{tab:q1-temperature} 和表~\ref{tab:q1-moisture}，结果保留四位小数。

\input{tables/q1_temperature}
\input{tables/q1_moisture}

药材表面先升温、失水，中心响应滞后。$1800\,\mathrm{s}$ 时中心与表面温差为 $@@TEMP_GAP@@\,{}^\circ\mathrm{C}$，二者仍低于烘房温度 $@@ENV_END@@\,{}^\circ\mathrm{C}$；表面含水率已降至 $@@SURFACE_C@@\,\mathrm{kg/kg}$，而中心四位小数下仍为 $2.5500$，表明内外热湿梯度仍然显著。

\Needspace{8\baselineskip}
\noindent\textbf{（2）传递特征}

图~\ref{fig:q1-overview} 以每秒、每隔 $0.1\,\mathrm{cm}$ 的未舍入结果给出完整时空场，并用规定时刻的径向切片核对其趋势；绘图未另作拟合平滑。

\begin{figure}[!htbp]
\centering
\includegraphics[width=\textwidth]{figures/q1_overview.pdf}
\caption{预热阶段温度与干基含水率的时空演化及径向切片}
\label{fig:q1-overview}
\end{figure}

图~\ref{fig:q1-overview} 显示温度扰动已传至中心，而含水率低值区仍集中于外层，说明热量向内传递快于水分向外扩散；且 $C$ 降低使 $D(C)$ 减小，外层干燥后补水进一步受限。$1800\,\mathrm{s}$ 时体积加权平均含水率为 $@@MEAN_C@@\,\mathrm{kg/kg}$，较初值下降 $@@DROP@@\%$，药材尚未均匀干燥。
