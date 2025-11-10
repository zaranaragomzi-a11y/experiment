import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# 페이지 설정
st.set_page_config(page_title="가상 산-염기 적정 실험실 (확장 버전)", layout="wide")

st.title("🔬 가상 산-염기 적정 실험실 (확장 버전)")
st.markdown("**가상의 적정 곡선을 시각화하고, 실제 실험 데이터와 비교할 수 있습니다.**")

# ------------------------------
# 1️⃣ 산 / 염기 선택
# ------------------------------
st.sidebar.header("1️⃣ 실험 조건 선택")

acid_type = st.sidebar.selectbox("산의 종류", ["강산", "약산"])
base_type = st.sidebar.selectbox("염기의 종류", ["강염기", "약염기"])

acid_eq = st.sidebar.selectbox("산의 당량수 (가수)", [1, 2, 3])
base_eq = st.sidebar.selectbox("염기의 당량수 (가수)", [1, 2, 3])

acid_conc = st.sidebar.number_input("산의 농도 (M)", 0.1, 2.0, 0.1, step=0.1)
base_conc = st.sidebar.number_input("염기의 농도 (M)", 0.1, 2.0, 0.1, step=0.1)
acid_vol = st.sidebar.number_input("산의 부피 (mL)", 10.0, 100.0, 25.0, step=5.0)

# 약산/약염기 Ka, Kb 값 입력
if acid_type == "약산":
    Ka = 10 ** (-st.sidebar.number_input("약산의 pKa", 3.0, 10.0, 5.0))
else:
    Ka = None
if base_type == "약염기":
    Kb = 10 ** (-st.sidebar.number_input("약염기의 pKb", 3.0, 10.0, 5.0))
else:
    Kb = None

# ------------------------------
# 2️⃣ 지시약 선택
# ------------------------------
st.sidebar.header("2️⃣ 지시약 선택")

indicators = {
    "메틸 오렌지": (3.1, 4.4),
    "메틸 레드": (4.4, 6.2),
    "브로모티몰 블루": (6.0, 7.6),
    "페놀프탈레인": (8.3, 10.0),
    "티몰 블루": (1.2, 2.8)
}
indicator_name = st.sidebar.selectbox("지시약 선택", list(indicators.keys()))
indicator_range = indicators[indicator_name]

# ------------------------------
# 3️⃣ CSV 데이터 업로드
# ------------------------------
st.sidebar.header("3️⃣ 실험 데이터 비교")
uploaded_file = st.sidebar.file_uploader("pH 데이터 CSV 업로드 (부피[mL], pH)", type=["csv"])

# ------------------------------
# 4️⃣ 계산 함수 정의
# ------------------------------
def calc_pH(Vb):
    Ca, Cb = acid_conc, base_conc
    Va = acid_vol
    nA = Ca * Va / 1000 * acid_eq
    nB = Cb * Vb / 1000 * base_eq

    # 강산-강염기
    if acid_type == "강산" and base_type == "강염기":
        if nB < nA:
            H = (nA - nB) / ((Va + Vb) / 1000)
            pH = -np.log10(H)
        elif nB == nA:
            pH = 7.0
        else:
            OH = (nB - nA) / ((Va + Vb) / 1000)
            pH = 14 + np.log10(OH)
        return pH

    # 약산-강염기
    elif acid_type == "약산" and base_type == "강염기":
        Ka_local = Ka
        if nB < nA:
            nHA = nA - nB
            nA_minus = nB
            pH = np.log10(nA_minus / nHA) + (-np.log10(Ka_local))
        elif nB == nA:
            pH = 7 + 0.5 * (14 + np.log10(Ka_local))
        else:
            OH = (nB - nA) / ((Va + Vb) / 1000)
            pH = 14 + np.log10(OH)
        return pH

    # 강산-약염기
    elif acid_type == "강산" and base_type == "약염기":
        Kb_local = Kb
        if nB > nA:
            OH = (nB - nA) / ((Va + Vb) / 1000)
            pH = 14 + np.log10(OH)
        elif nB == nA:
            pH = 7 - 0.5 * (14 + np.log10(Kb_local))
        else:
            nB_remaining = nA - nB
            nBH_plus = nB
            pH = 14 - (-np.log10(Kb_local) + np.log10(nB_remaining / nBH_plus))
        return pH

    else:
        return 7.0

# ------------------------------
# 5️⃣ 적정 곡선 계산
# ------------------------------
Vb_values = np.linspace(0, 2 * acid_vol, 300)
pH_values = [calc_pH(v) for v in Vb_values]
diffs = np.gradient(pH_values)
eq_index = np.argmax(diffs)
eq_vol = Vb_values[eq_index]
eq_pH = pH_values[eq_index]

# ------------------------------
# 6️⃣ 그래프 시각화
# ------------------------------
st.header("📈 적정 곡선 시각화")

fig = go.Figure()

# (1) 시뮬레이션 곡선
fig.add_trace(go.Scatter(
    x=Vb_values, y=pH_values, mode='lines', name='시뮬레이션 곡선',
    line=dict(color='blue')
))

# (2) 변곡점 표시
fig.add_vline(x=eq_vol, line=dict(color='red', dash='dot'))
fig.add_annotation(x=eq_vol, y=eq_pH, text=f"중화점 ≈ {eq_vol:.1f} mL, pH={eq_pH:.2f}",
                   showarrow=True, arrowhead=2, arrowcolor='red')

# (3) 지시약 변색 범위 시각화
fig.add_vrect(
    x0=min(Vb_values), x1=max(Vb_values),
    y0=indicator_range[0], y1=indicator_range[1],
    fillcolor="green", opacity=0.15,
    annotation_text=f"{indicator_name} ({indicator_range[0]}~{indicator_range[1]})",
    annotation_position="top left"
)

# (4) 업로드된 실제 데이터 표시
if uploaded_file is not None:
    df_exp = pd.read_csv(uploaded_file)
    fig.add_trace(go.Scatter(
        x=df_exp.iloc[:, 0], y=df_exp.iloc[:, 1],
        mode='markers+lines', name='실험 데이터', marker=dict(color='orange', size=6)
    ))

# 그래프 설정
fig.update_layout(
    xaxis_title="염기 부피 (mL)",
    yaxis_title="pH",
    template="plotly_white",
    width=900, height=550,
    legend=dict(x=0.02, y=0.98)
)

st.plotly_chart(fig, use_container_width=True)

# ------------------------------
# 7️⃣ 결과 분석
# ------------------------------
st.markdown(f"**⚗️ 중화점:** 약 {eq_vol:.2f} mL, pH = {eq_pH:.2f}")

if indicator_range[0] <= eq_pH <= indicator_range[1]:
    st.success(f"✅ {indicator_name}는 적합한 지시약입니다.")
else:
    st.error(f"❌ {indicator_name}는 적합하지 않습니다. 변곡점 pH = {eq_pH:.2f}")

if acid_type == "약산":
    half_eq_vol = eq_vol / 2
    half_pH = calc_pH(half_eq_vol)
    st.info(f"📊 pKa ≈ {half_pH:.2f}")

# ------------------------------
# 8️⃣ 조합 비교 모드
# ------------------------------
st.header("🔍 여러 조건 비교 모드")

if "saved_curves" not in st.session_state:
    st.session_state.saved_curves = []

if st.button("현재 설정 곡선 저장"):
    st.session_state.saved_curves.append({
        "acid": acid_type, "base": base_type,
        "curve_x": Vb_values, "curve_y": pH_values
    })
    st.success("현재 조합의 곡선이 저장되었습니다.")

if len(st.session_state.saved_curves) > 0:
    fig_compare = go.Figure()
    for i, curve in enumerate(st.session_state.saved_curves):
        fig_compare.add_trace(go.Scatter(
            x=curve["curve_x"], y=curve["curve_y"],
            mode='lines', name=f"{curve['acid']}-{curve['base']} #{i+1}"
        ))
    fig_compare.update_layout(
        title="조합 비교 곡선",
        xaxis_title="염기 부피 (mL)",
        yaxis_title="pH",
        template="plotly_white"
    )
    st.plotly_chart(fig_compare, use_container_width=True)
