import streamlit as st

st.set_page_config(layout="wide")

st.components.v1.iframe(
    src="http://127.0.0.1:3838",
    width=1400,
    height=900,
    scrolling=True
)
