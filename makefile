# 定义编译器
CXX = g++

# 定义编译选项
CXXFLAGS = -std=c++11 -Wall -O2

# 定义目标文件
TARGET = bin/FC-Virus

# 定义源文件和头文件的路径
SRCS = code/main.cpp code/exl_1.cpp code/ex_r.cpp code/GeneralSet.cpp code/kmer.cpp code/loadreads.cpp
HEADERS = code/ex_l.h code/ex_r.h code/GeneralSet.h code/kmer.h code/loadreads.h

# 定义目标规则
OBJS = $(SRCS:code/%.cpp=code/%.o)

# 链接目标文件到 bin 目录
$(TARGET): $(OBJS)
	@mkdir -p bin  # 如果 bin 目录不存在，则创建
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# 编译每个.cpp文件到 code 目录中的 .o 文件
code/%.o: code/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 清除生成的目标文件和可执行文件
clean:
	rm -f code/*.o $(TARGET)

.PHONY: clean
