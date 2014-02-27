{-# OPTIONS -fglasgow-exts #-}
module Tiff.LibTIFF -- (open, setField, tIFFTAG_IMAGEWIDTH .....
where

-- typedef	uint32 ttag_t;		/* directory tag */
-- typedef	uint16 tdir_t;		/* directory index */
-- typedef	uint16 tsample_t;	/* sample number */
-- typedef	uint32 tstrip_t;	/* strip number */
-- typedef uint32 ttile_t;		/* tile number */
-- typedef	int32 tsize_t;		/* i/o size in bytes */
-- typedef	void* tdata_t;		/* image data ref */
-- typedef	uint32 toff_t;		/* file offset */

import Foreign.Ptr
import Foreign.C
import Foreign.C.Types
import Foreign.Marshal.Array

type TIFF = Ptr ()

foreign import ccall unsafe "TIFFOpen" tIFFOpen :: CString -> CString -> IO TIFF

-- 	tsize_t  TIFFWriteEncodedStrip(TIFF*  tif,  tstrip_t  strip,
-- 	tdata_t buf, tsize_t size)
foreign import ccall unsafe "TIFFWriteEncodedStrip" tIFFWriteEncodedStrip :: TIFF -> CUInt -> (Ptr ()) -> CUInt -> IO CUInt

foreign import ccall unsafe "TIFFClose" tIFFClose :: TIFF -> IO ()

-- ==============================================================
--      high level functions  -- wrappers
open :: String -> String -> IO TIFF
open filename mode = withCString filename $ \f -> withCString mode $ \m -> tIFFOpen f m
-- ==============================================================


type Tag = CInt

-- 	TIFFTAG_IMAGEWIDTH              1      uint32            
tIFFTAG_IMAGEWIDTH :: TIFF -> CUInt -> IO ()
-- curry
tIFFTAG_IMAGEWIDTH image value = setUint32 image 256 value

foreign import ccall unsafe "TIFFSetField" setUint32 :: TIFF -> Tag -> CUInt -> IO ()
foreign import ccall unsafe "TIFFSetField" setUint16 :: TIFF -> Tag -> CUShort -> IO ()
foreign import ccall unsafe "TIFFSetField" setDouble :: TIFF -> Tag -> CDouble -> IO ()


type TagF a = TIFF -> a -> IO ()

setField :: TIFF -> (TagF a) -> a -> IO ()
setField image f value = f image value

-- ================

-- setFields ??  

-- Template Haskell??

tIFFTAG_IMAGELENGTH, tIFFTAG_ROWSPERSTRIP :: TIFF -> CUInt -> IO ()
tIFFTAG_IMAGELENGTH i v = setUint32 i 257 v

-- |
-- | StripsPerImage = floor ((ImageLength + RowsPerStrip - 1) / RowsPerStrip).
-- | StripsPerImage is not a field. It is merely a value that a TIFF reader will want to compute because it specifies the number of StripOffsets and StripByteCounts for the image. 
tIFFTAG_ROWSPERSTRIP i v = setUint32 i 278 v

tIFFTAG_BITSPERSAMPLE, tIFFTAG_SAMPLESPERPIXEL, tIFFTAG_COMPRESSION, tIFFTAG_PHOTOMETRIC, tIFFTAG_FILLORDER, tIFFTAG_PLANARCONFIG, tIFFTAG_RESOLUTIONUNIT :: TIFF -> CUShort -> IO ()
tIFFTAG_BITSPERSAMPLE i v = setUint16 i 258 v
tIFFTAG_SAMPLESPERPIXEL i v = setUint16 i 277 v
tIFFTAG_COMPRESSION i v = setUint16 i 259 v

cOMPRESSION_CCITTFAX4, cOMPRESSION_DEFLATE, cOMPRESSION_NONE, pHOTOMETRIC_MINISWHITE, pHOTOMETRIC_MINISBLACK, eXTRASAMPLE_UNASSALPHA, fILLORDER_MSB2LSB, pLANARCONFIG_CONFIG, rESUNIT_INCH, pHOTOMETRIC_RGB :: CUShort
cOMPRESSION_CCITTFAX4  = 4
cOMPRESSION_DEFLATE = 32946	-- Deflate compression
cOMPRESSION_NONE = 1 -- dump mode

tIFFTAG_PHOTOMETRIC i v = setUint16 i 262 v
pHOTOMETRIC_MINISWHITE = 0
pHOTOMETRIC_MINISBLACK = 1 -- min value is black
pHOTOMETRIC_RGB = 2

tIFFTAG_FILLORDER i v = setUint16 i 266 v
fILLORDER_MSB2LSB = 1
tIFFTAG_PLANARCONFIG i v = setUint16 i 284 v
pLANARCONFIG_CONFIG = 1
tIFFTAG_RESOLUTIONUNIT  i v = setUint16 i 296 v
rESUNIT_INCH = 2

tIFFTAG_XRESOLUTION, tIFFTAG_YRESOLUTION :: TIFF -> CDouble -> IO ()
tIFFTAG_XRESOLUTION i v = setDouble i 282 v
tIFFTAG_YRESOLUTION i v = setDouble i 283 v


foreign import ccall unsafe "TIFFSetField" setCUShortN' :: TIFF -> Tag -> Int -> (Ptr CUShort) -> IO ()
setCUShortN :: TIFF -> Tag -> [CUShort] -> IO ()
setCUShortN i t xs' = withArray xs' $ \xs -> setCUShortN' i t (length xs') xs 

tIFFTAG_EXTRASAMPLES :: TIFF -> [CUShort] -> IO ()
tIFFTAG_EXTRASAMPLES i v = setCUShortN i 338 v

-- | Unassociated alpha channels can be used to encode a number of independent masks, for example.
eXTRASAMPLE_UNASSALPHA = 2

tIFFTAG_XPOSITION, tIFFTAG_YPOSITION :: TIFF -> CDouble -> IO ()
tIFFTAG_XPOSITION i v = setDouble i 286 v
-- | The Y offset in ResolutionUnits of the top of the image, with respect to the top of the page. In the TIFF coordinate scheme, the positive Y direction is down, so that YPosition is always positive.
tIFFTAG_YPOSITION i v = setDouble i 287 v
