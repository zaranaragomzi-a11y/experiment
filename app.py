# fix syntax error for st.info string join
# replaced unterminated string with proper join

# snippet replacement
# in section for st.info of eq_points
# replace st.info("\n".join([...]))
# ensure proper escaping

# (We'll just output corrected version)

# === correction ===
if len(eq_points) > 0:
    st.subheader("중화점(등가점) 정보")
    rows = []
    for i, Veq in enumerate(eq_points):
        e_pH = eq_pH[i] if i < len(eq_pH) else np.nan
        ok = (ind_low <= e_pH <= ind_high)
        verdict = "✅ 적합" if ok else "⚠️ 부적합"
        msg = f"중화점의 pH는 {e_pH:.2f}이므로 {ind_name}은(는) {'적합합니다.' if ok else '적합하지 않습니다.'}"
        rows.append((i+1, Veq, e_pH, verdict, msg))
    st.table({
        "EQ#": [r[0] for r in rows],
        "V_eq (mL)": [f"{r[1]:.1f}" for r in rows],
        "pH_eq": [f"{r[2]:.2f}" for r in rows],
        "판정": [r[3] for r in rows]
    })
    st.info("\n".join([r[4] for r in rows]))
else:
    st.info("등가점을 계산할 수 없는 설정입니다. 농도/부피/종류를 확인하세요.")
