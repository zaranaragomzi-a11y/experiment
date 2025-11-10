import streamlit as st
import numpy as np
import plotly.graph_objects as go

st.title("💧 가상 산·염기 적정 실험실")

st.sidebar.header("🧪 실험 조건 입력")
acid = st.sidebar.selectbox("분석 물질(산)", ["HCl", "CH3COOH", "H2SO4"])
base = st.sidebar.selectbox("적정 용액(염기)", ["NaOH", "Ca(OH)2", "NH3"])
Ca = st.sidebar.number_input("산 농도 (M)", 0.01, 1.0, 0.1)
Va = st.sidebar.number_input("산 부피 (mL)", 1.0, 100.0, 25.0)
Cb = st.sidebar.number_input("염기 농도 (M)", 0.01, 1.0, 0.1)

st.write(f"선택된 조건: {acid}({Ca}M, {Va}mL) + {base}({Cb}M)")

# 간단한 pH 계산 예시 (강산-강염기 가정)
Vb = np.linspace(0, 50, 200)
nH = Ca * Va / 1000
pH = []
for V in Vb:
    nOH = Cb * V / 1000
    if nH > nOH:
        h = (nH - nOH) / ((Va + V) / 1000)
        pH.append(-np.log10(h))
    elif nH < nOH:
        oh = (nOH - nH) / ((Va + V) / 1000)
        pH.append(14 + np.log10(oh))
    else:
        pH.append(7)

# 그래프 표시
fig = go.Figure()
fig.add_trace(go.Scatter(x=Vb, y=pH, mode="lines", name="적정 곡선"))
fig.update_layout(xaxis_title="부피 (mL)", yaxis_title="pH", title="적정 곡선")
st.plotly_chart(fig, use_container_width=True)
