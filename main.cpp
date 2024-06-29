#include <iostream>
#include <conio.h>
#include <stdlib.h>
//#include <string>

int main(void) {
	int screen[30][30];
	std::string printString = "jkl;jkljkl;jkl;jkl;jkl;jkl;";

	for (int i = 0; i < 30; i++) {
		for (int j = 0; j < 30; j++) {
			screen[i][j] = 0;
		}
	}

	while (1 == 1) {
		for (int i = 0; i < 30; i++) {
			for (int j = 0; j < 30; j++) {
				std::cout << screen[i][j] << " ";
			}
			std::cout << std::endl;
		}
		system("cls");

		for (int i = 0; i < 30; i++) {
			for (int j = 0; j < 30; j++) {
				std::cout << "1" << " ";
			}
			std::cout << std::endl;
		}
		system("cls");
	}

	return 0;
}