@echo off

rem create virtual environment using Python 3.9
rem python -m venv venv --python=python3.9

rem activate virtual environment
call venv\Scripts\activate

pip freeze > requirements.txt

echo Completed!
pause

