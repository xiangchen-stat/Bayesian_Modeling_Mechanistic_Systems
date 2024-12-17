#include <stdio.h>

void __stack_chk_fail() {
    printf("Stack check failed!\n");
}
