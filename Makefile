	
CXX = g++
CXXFLAGS = -O3 -Wall -Wshadow -Werror -std=c++11 -pedantic

INCLUDES =
LDFLAGS =
LIBS =

TARGET = lbm
OBJS = $(TARGET).o GrayScaleImage.o lodepng.o

all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

$(TARGET).o: $(TARGET).cpp Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(TARGET).cpp

GrayScaleImage.o: imageClass/GrayScaleImage.cpp imageClass/GrayScaleImage.h Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<
	
lodepng.o: imageClass/lodepng.cpp imageClass/lodepng.h Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<

clean:
	@$(RM) -rf *.o $(TARGET)
