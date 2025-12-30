        im = protein_map_series(:,:,step); % double precision
        t = Tiff(strcat('D:\Kaustav\CODES SET\Membrane_Sheet_Simulation\dynamic\06\results_protein_map\protein_map', sprintf('%05d',step),'.tif'), 'w');
        
        tagstruct.ImageLength = size(im, 1);
        tagstruct.ImageWidth = size(im, 2);
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 64;  % 64-bit float
        tagstruct.SamplesPerPixel = 1;
        tagstruct.RowsPerStrip = 16;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.Compression = Tiff.Compression.None;
        
        t.setTag(tagstruct);
        t.write(im);  % No casting needed, already double
        t.close();