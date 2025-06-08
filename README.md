# simpleFedback
Simple Feedback LV2 plugin

compile with
gcc -fvisibility=hidden -fPIC -Wl,-Bstatic -Wl,-Bdynamic -Wl,--as-needed -shared -pthread `pkg-config --cflags lv2` -lm `pkg-config --libs lv2` simpleFeedback.c -o simpleFeedback.so
