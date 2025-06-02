@echo off

rem activate virtual environment
call venv\Scripts\activate

pip install -r requirements.txt

echo Completed!
pause
