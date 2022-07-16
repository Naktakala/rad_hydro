ffmpeg -f image2 -framerate 5 \
-i Test3a_t%003d.png \
-pix_fmt rgb8 \
-vf scale=1280:-1 \
Y2Video.mp4

ffmpeg -i Y2Video.mp4 -vf "scale=1280:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" Y2Gif.gif

#-i Test3a_t002.png \
#-i Test3a_t003.png \
#-i Test3a_t004.png \