@echo off
echo Loading packages
call venv\Scripts\activate
echo Launching Graphical Interface
streamlit run dash.py
pause