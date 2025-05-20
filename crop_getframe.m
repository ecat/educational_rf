function [getframe_struct] =  crop_getframe(getframe_struct, xlims, ylims)

getframe_struct.cdata = getframe_struct.cdata(xlims, ylims, :);


