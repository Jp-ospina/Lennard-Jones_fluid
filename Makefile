CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic
CPPFLAGS ?=
LDFLAGS ?=

TARGET := fluido.exe
SOURCES := src/main.cpp src/main_sin_modularizar.cpp
HEADERS := include/fluido.h

.PHONY: all run clean

ifeq ($(OS),Windows_NT)
REMOVE_FILE := cmd.exe /C del /Q
else
REMOVE_FILE := rm -f
endif

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude $(SOURCES) $(LDFLAGS) -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	$(REMOVE_FILE) $(TARGET)
