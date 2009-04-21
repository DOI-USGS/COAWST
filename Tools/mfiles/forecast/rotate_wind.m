function [windu_rot, windv_rot]=rotate_wind(ugridcomb, vgridcomb, roms_angle)

windu_rot=ugridcomb.*cos(roms_angle(1,1))+vgridcomb.*sin(roms_angle(1,1)); 
windv_rot=vgridcomb.*cos(roms_angle(1,1))-ugridcomb.*sin(roms_angle(1,1)); 