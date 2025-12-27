# gcc -fvisibility=hidden -fPIC -Wl,-Bstatic -Wl,-Bdynamic -Wl,--as-needed -shared -pthread `pkg-config --cflags lv2` -lm `pkg-config --libs lv2` simpleFeedback.c -o simpleFeedback.so
CFLAGS += -DDEBUG
all: simpleFeedback.c
	gcc -fvisibility=hidden -fPIC -Wl,-Bstatic -Wl,-Bdynamic -Wl,--as-needed -shared -pthread `pkg-config --cflags lv2` -lm `pkg-config --libs lv2` simpleFeedback.c -o simpleFeedback.so -lfftw3f
	mv simpleFeedback.so /home/$(USER)/.lv2/simpleFeedback/simpleFeedback.so