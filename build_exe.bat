@echo off
chcp 65001 >nul
echo ============================================
echo   基因型分析 - 打包为 Windows exe
echo ============================================
echo.

cd /d "%~dp0"

REM 检查 Python
python --version >nul 2>&1
if errorlevel 1 (
    echo [错误] 未找到 Python，请先安装 Python 并加入 PATH。
    pause
    exit /b 1
)

REM 安装 PyInstaller（若未安装）
pip show pyinstaller >nul 2>&1
if errorlevel 1 (
    echo 正在安装 PyInstaller ...
    pip install pyinstaller
    if errorlevel 1 (
        echo [错误] 安装 PyInstaller 失败。
        pause
        exit /b 1
    )
)

REM 确保依赖已安装
echo 检查依赖 ...
pip install -r requirements.txt -q
pip install parasail -q

REM 清理旧构建
if exist "dist\GeneType.exe" del "dist\GeneType.exe"
if exist "build" rmdir /s /q build

REM 使用 spec 打包（单文件、无控制台）
echo.
echo 正在打包，请稍候 ...
pyinstaller --noconfirm GeneType.spec

if errorlevel 1 (
    echo.
    echo [错误] 打包失败。
    pause
    exit /b 1
)

echo.
echo ============================================
echo   打包完成
echo ============================================
echo   可执行文件: dist\GeneType.exe
echo   可直接复制 GeneType.exe 到任意 Windows 10/11 电脑运行。
echo   无需安装 Python。
echo ============================================
echo.
explorer dist
pause
