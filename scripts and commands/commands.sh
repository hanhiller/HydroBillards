
# -------------- Movie making and minor editing (using ffmpeg) commands --------------
# Movie creation with image starting from 'frame_0.png'
ffmpeg -framerate 25 -i frame_%d.png -c:v libx265 -pix_fmt yuv420p movie.mp4

# Movie creation with image starting from 'frame_2291.png'
ffmpeg -framerate 15 -start_number 2291 -i frame_%d.png -c:v libx265 -pix_fmt yuv420p movie.mp4

# Movie creation with image starting from 'frame_2291.png' and upto certain frames
ffmpeg -framerate 2 -start_number 2291 -i frame_%d.png -vframes 100 -c:v libx265 -pix_fmt yuv420p movie.mp4

# Cut movie
ffmpeg -i movie.mp4 -vcodec libx265 -crf 34 com_movie.mp4

# find video resolution
ffprobe -v error -select_streams v:0 -show_entries stream=width,height -of default=nw=1 movie_T60_F35.mp4

# Rescale a video
ffmpeg -i movie_T60_F35.mp4 -vf scale=1920:800 movie_T60_F35_rescaled.mp4

# Double the speed of the video
ffmpeg -i movie_T60_F35.mp4 -filter:v "setpts=0.5*PTS" movie_T60_2x.mp4



# ---------------------------- Simulation specific commands --------------------------

# running simulation command
python3 main.py simulate --sim-type=full_circle --temperature=90 --time-steps=5000 --output-path=.