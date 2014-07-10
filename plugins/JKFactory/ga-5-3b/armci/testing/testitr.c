#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#include "iterator.h"

int main(int argc, char **argv)
{
  stride_info_t sinfo, dinfo;
  int a[10][10], b[11][11];
  int asr[1] = {10 * sizeof(int)};
  int bsr[1] = {11 * sizeof(int)};
  int count[2] = {5 * sizeof(int), 5};
  int i, j;

  for (i = 0; i < 10; i++) {
    for (j = 0; j < 10; j++) {
      a[i][j] = i * 10 + j;
    }
  }
  for (i = 0; i < 11; i++) {
    for (j = 0; j < 11; j++) {
      b[i][j] = 0;
    }
  }

  armci_stride_info_init(&sinfo, &a[2][3], 1, asr, count);
  armci_stride_info_init(&dinfo, &b[3][4], 1, bsr, count);

  assert(armci_stride_info_size(&sinfo) == 5);
  assert(armci_stride_info_size(&dinfo) == 5);
  assert(armci_stride_info_pos(&sinfo) == 0);
  assert(armci_stride_info_pos(&dinfo) == 0);

  while (armci_stride_info_has_more(&sinfo)) {
    int bytes;
    char *ap, *bp;
    assert(armci_stride_info_has_more(&dinfo));

    bytes = armci_stride_info_seg_size(&sinfo);
    assert(bytes == armci_stride_info_seg_size(&dinfo));

    ap = armci_stride_info_seg_ptr(&sinfo);
    bp = armci_stride_info_seg_ptr(&dinfo);

    memcpy(bp, ap, bytes);

    armci_stride_info_next(&sinfo);
    armci_stride_info_next(&dinfo);
  }
  armci_stride_info_destroy(&sinfo);
  armci_stride_info_destroy(&dinfo);

#if 0
  for (i = 0; i < 10; i++) {
    for (j = 0; j < 10; j++) {
      printf("%3d  ", a[i][j]);
    }
    printf("\n");
  }

  for (i = 0; i < 11; i++) {
    for (j = 0; j < 11; j++) {
      printf("%3d  ", b[i][j]);
    }
    printf("\n");
  }
#endif

  for (i = 2; i < 2 + 5; i++) {
    for (j = 3; j < 3 + 5; j++) {
      if (a[i][j] != b[i+1][j+1]) {
        printf("a[%d][%d]=%d b[%d][%d]=%d\n", i, j, a[i][j], i, j, b[i][j]);
        printf("Test Failed\n");
        return 0;
      }
    }
  }
  printf("Success\n");
  return 0;
}

