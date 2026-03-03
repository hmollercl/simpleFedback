
#all: simpleFeedback.c
#	gcc -fvisibility=hidden -fPIC -Wl,-Bstatic -Wl,-Bdynamic -Wl,--as-needed -shared -pthread `pkg-config --cflags lv2` -lm `pkg-config --libs lv2` simpleFeedback.c -o simpleFeedback.so -lfftw3f
#	mv simpleFeedback.so /home/$(USER)/.lv2/simpleFeedback/simpleFeedback.so

CC      ?= gcc
CFLAGS  += -fvisibility=hidden -fPIC -pthread -DDEBUG `pkg-config --cflags lv2`
LDFLAGS += -shared -Wl,-Bstatic -Wl,-Bdynamic -Wl,--as-needed
LDLIBS  += -lm -lfftw3f `pkg-config --libs lv2`

TARGET  = simpleFeedback.so
SRC     = simpleFeedback.c
INSTALL_DIR = /home/$(USER)/.lv2/simpleFeedback

all: $(TARGET)
	mkdir -p $(INSTALL_DIR)
	mv $(TARGET) $(INSTALL_DIR)/$(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(TARGET)