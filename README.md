# simpleFedback
Simple Feedback LV2 plugin
currently it's only is a small mono delay
TODO:
detect frequency and add filter for it.

to build:
```gcc -fvisibility=hidden -fPIC -Wl,-Bstatic -Wl,-Bdynamic -Wl,--as-needed -shared -pthread `pkg-config --cflags lv2` -lm `pkg-config --libs lv2` channelSwitch.c -o channelSwitch.so```

needed in ubuntu:
```sudo apt install build-essential pkg-config lv2-dev```
