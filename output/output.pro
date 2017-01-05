TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH=/home/huheng/source/ffmpeg
LIBS=-lavcodec -lavfilter -lavutil -lswresample \
-lavformat -lswscale  -lavdevice -lpostproc -lpthread




SOURCES += main.c
