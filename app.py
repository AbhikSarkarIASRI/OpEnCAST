import streamlit as st

st.set_page_config(layout="wide")
st.title("OpEnCAST")

st.markdown(
    '<iframe src="http://localhost:3838" width="100%" height="900"></iframe>',
    unsafe_allow_html=True
)
