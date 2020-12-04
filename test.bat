echo 正在清理VS2015工程中不需要的文件
echo 开始清理请稍等......
@echo off
 
rem 清理
attrib -s -h -r *opendb
del /s /q  *.opendb
del /s /q  *.db
 
del /s /q .\ipch\*.*
rd  /s /q .\ipch\
 
del /s /q .\Debug\*.*
rd  /s /q .\Debug\
 
del /s /q .\Release\*.*
rd  /s /q .\Release\
 
 
for %%p in ("%cd%") do set folder=%%~nxp
echo %folder%
del /s /q .\%folder%\Debug\*.*
rd  /s /q .\%folder%\Debug\
del /s /q .\%folder%\Release\*.*
rd  /s /q .\%folder%\Release\