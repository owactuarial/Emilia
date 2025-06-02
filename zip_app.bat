@echo off

rem activate virtual environment
call venv\Scripts\activate

python zip_app.py

echo Completed!
pause