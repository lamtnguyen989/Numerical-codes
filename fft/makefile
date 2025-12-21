SRC_DIR=src
BIN_DIR=bin

basic: $(SRC_DIR)/fft.cc
	g++ $< -O3 -flto -o $(BIN_DIR)/fft

tensor: $(SRC_DIR)/tensor.h $(SRC_DIR)/fft_tensor.cc
	g++ -O3 -flto $< -o $(BIN_DIR)/fft_tensor