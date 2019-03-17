% load the LIC
lic = readNPY('./lic.npy');
% 250 um
I250 = fitsread('carinaData/blastData/carinaneb_good_250_p10_good_C_gls_map_cal.fits', 'image', 1);
Q250 = fitsread('carinaData/blastData/carinaneb_good_250_p10_good_C_gls_map_cal.fits', 'image', 2);
U250 = fitsread('carinaData/blastData/carinaneb_good_250_p10_good_C_gls_map_cal.fits', 'image', 3);

% 350 um
I300 = fitsread('carinaData/blastData/carinaneb_good_350_p10_good_C_gls_map_cal.fits', 'image', 1);
Q300 = fitsread('carinaData/blastData/carinaneb_good_350_p10_good_C_gls_map_cal.fits', 'image', 2);
U300 = fitsread('carinaData/blastData/carinaneb_good_350_p10_good_C_gls_map_cal.fits', 'image', 3);

% 500 um
I500 = fitsread('carinaData/blastData/carinaneb_good_500_p10_good_C_gls_map_cal.fits', 'image', 1);
Q500 = fitsread('carinaData/blastData/carinaneb_good_500_p10_good_C_gls_map_cal.fits', 'image', 2);
U500 = fitsread('carinaData/blastData/carinaneb_good_500_p10_good_C_gls_map_cal.fits', 'image', 3);

figure;
subplot(3,1,1);imagesc(I250);
subplot(3,1,2);imagesc(I300);
subplot(3,1,3);imagesc(I500);
figure;
subplot(3,1,1); edge250 = edge(I250,'Canny'); imagesc(edge250); colormap(copper);
hold on;
subplot(3,1,2); edge300 = edge(I300,'Canny'); imagesc(edge300); colormap(copper)
hold on;
subplot(3,1,3); edge500 = edge(I500,'Canny'); imagesc(edge500); colormap(copper)
hold off;
