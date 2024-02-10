#include <raylib.h>

static int do_main() {
  InitWindow(600, 600, "a");
  SetTargetFPS(90);
  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground(BLACK);
    EndDrawing();
  }
  CloseWindow();
  return 0;
}

#if defined(_WIN32)

int wWinMain(void **_0, void **_1, void **_2, int _3) { return do_main(); }

#else

int main() { return do_main(); }

#endif
