gif-h
=====

This one-header library offers a simple, very limited way to create animated GIFs directly in code.
There's not much effort to optimize colors to improve compression, so resulting files
are often quite large. But it's still handy to be able to spit out GIFs on command.

Only RGBA8 is currently supported as an input format; the alpha is ignored.

USAGE:
Create a GifWriter struct. Pass it to GifBegin() to initialize and write the header.
Pass subsequent frames to GifWriteFrame().
Finally, call GifEnd() to close the file handle and free memory.
