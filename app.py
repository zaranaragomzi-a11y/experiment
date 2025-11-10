import streamlit as st
import numpy as np
import pandas as pd
from PIL import Image
import requests
from io import BytesIO
import plotly.graph_objects as go
from time import sleep

# ------------------------------
# 페이지 설정
# ------------------------------
st.set_page_config(page_title="교육용 산-염기 적정 실험실", layout="wide")
st.title("🏫 가상 산-염기 적정 실험실 (교육용)")

# ------------------------------
# GitHub raw URL로 이미지 불러오기
# ------------------------------
buret_url = "https://raw.githubusercontent.com/zaranaragomzi-a11y/experiment/main/images/buret.png"
flask_url = "https://raw.githubusercontent.com/zaranaragomzi-a11y/experiment/main/images/flask.png"

buret_img = Image.open(BytesIO(requests.get(buret_url).content))
flask_img = Image.open(BytesIO(requests.get(flask_url).content))

# ------------------------------
# 실험 조건
# ------------------------------
st.sidebar.header("실험 조건")
acid_type = st.sidebar.selectbox("산 종류", ["강산", "약산"])
base_type = st.sidebar.selectbox("염기 종류", ["강염기", "약염기"])
acid_eq = st.sidebar.selectbox("산 당량수", [1,2,3])
base_eq = st.sidebar.selectbox("염기 당량수", [1,2,3])
acid_conc = st.sidebar.number_input("산 농도 (M)", 0.1, 2.0, 0.1)
acid_vol = st.sidebar.number_input("산 부피 (mL)", 10.0, 100.0, 25.0)
base_conc = st.sidebar.number_input("염기 농도 (M)", 0.1, 2.0, 0.1)
base_vol = st.sidebar.slider("적정 용액 최대 부피 (mL)", 0.0, 2*acid_vol, 50.0)

Ka = 10**(-st.sidebar.number_input("약산 pKa",3.0,10.0,5.0)) if acid_type=="약산" else None
Kb = 10**(-st.sidebar.number_input("약염기 pKb",3.0,10.0,5.0)) if base_type=="약염기" else None

# ------------------------------
# 지시약 선택
# ------------------------------
indicators = {
    "메틸 오렌지": (3.1, 4.4,"#ff9900"),
    "메틸 레드": (4.4,6.2,"#ff0000"),
    "브로모티몰 블루": (6.0,7.6,"#00ff00"),
    "페놀프탈레인": (8.3,10.0,"#9900ff"),
    "티몰 블루": (1.2,2.8,"#0000ff")
}
indicator_name = st.sidebar.selectbox("지시약 선택", list(indicators.keys()))
indicator_range, indicator_color = indicators[indicator_name][:2], indicators[indicator_name][2]

# ------------------------------
# CSV 업로드
# ------------------------------
uploaded_file = st.sidebar.file_uploader("실험 데이터 CSV", type=["csv"])

# ------------------------------
# pH 계산 함수
# ------------------------------
def calc_pH(Vb):
    nA = acid_conc*acid_vol/1000*acid_eq
    nB = base_conc*Vb/1000*base_eq
    if acid_type=="강산" and base_type=="강염기":
        if nB<nA: return -np.log10((nA-nB)/((acid_vol+Vb)/1000))
        elif nB==nA: return 7
        else: return 14+np.log10((nB-nA)/((acid_vol+Vb)/1000))
    elif acid_type=="약산" and base_type=="강염기":
        if nB<nA: return np.log10(nB/(nA-nB))-np.log10(Ka)
        elif nB==nA: return 7+0.5*(14+np.log10(Ka))
        else: return 14+np.log10((nB-nA)/((acid_vol+Vb)/1000))
    elif acid_type=="강산" and base_type=="약염기":
        if nB>nA: return 14+np.log10((nB-nA)/((acid_vol+Vb)/1000))
        elif nB==nA: return 7-0.5*(14+np.log10(Kb))
        else: return 14-(-np.log10(Kb)+np.log10(nA-nB)/nB)
    else: return 7

# ------------------------------
# 적정 시뮬레이션
# ------------------------------
Vb_values = np.linspace(0, base_vol, 50)
pH_values = [calc_pH(v) for v in Vb_values]
diffs = np.gradient(pH_values)
eq_index = np.argmax(diffs)
eq_vol, eq_pH = Vb_values[eq_index], pH_values[eq_index]

# 빈 Plotly 그래프 생성
fig = go.Figure()
fig.add_trace(go.Scatter(x=[], y=[], mode='lines+markers', name='적정 곡선', line=dict(color='blue')))
fig.update_layout(xaxis_title="염기 부피 (mL)", yaxis_title="pH", template="plotly_white", width=900, height=500)
plot_area = st.empty()

# 플라스크, 뷰렛 열
flask_col, buret_col = st.columns([1,0.2])
flask_disp = flask_col.empty()
buret_disp = buret_col.empty()
buret_disp.image(buret_img, width=50)

# ------------------------------
# 방울 단위 애니메이션
# ------------------------------
x_vals, y_vals = [], []
for V, pH in zip(Vb_values, pH_values):
    # Plotly 그래프 업데이트
    x_vals.append(V)
    y_vals.append(pH)
    fig.data[0].x = x_vals
    fig.data[0].y = y_vals
    plot_area.plotly_chart(fig, use_container_width=True)
    
    # 플라스크 색상 결정
    if indicator_range[0]<=pH<=indicator_range[1]:
        color = indicator_color
    elif pH<indicator_range[0]:
        color = "#ff6666" # 산성
    else:
        color = "#6666ff" # 염기성
    flask_disp.image(flask_img)
    flask_disp.markdown(f"<div style='width:100%;height:150px;background-color:{color};margin-top:-140px'></div>", unsafe_allow_html=True)
    
    sleep(0.1)

# ------------------------------
# 결과 표시
# ------------------------------
st.markdown(f"**중화점:** {eq_vol:.2f} mL, pH = {eq_pH:.2f}")
if indicator_range[0]<=eq_pH<=indicator_range[1]:
    st.success(f"{indicator_name} 적합")
else:
    st.error(f"{indicator_name} 부적합")

if acid_type=="약산":
    half_pH = calc_pH(eq_vol/2)
    st.info(f"pKa ≈ {half_pH:.2f}")

# ------------------------------
# CSV 데이터 표시
# ------------------------------
if uploaded_file:
    df_exp = pd.read_csv(uploaded_file)
    fig.add_trace(go.Scatter(x=df_exp.iloc[:,0], y=df_exp.iloc[:,1], mode='markers+lines', name='실험 데이터',
                             marker=dict(color='orange', size=6)))
    plot_area.plotly_chart(fig, use_container_width=True)
