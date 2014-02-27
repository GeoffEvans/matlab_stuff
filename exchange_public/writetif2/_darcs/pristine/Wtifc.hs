-- in order to run add C:\GnuWin32\bin to %PATH%
module Wtifc where

import Matlab.Mex
import qualified Matlab.Mex as M (map)
import Tiff.LibTIFF

-- XXX
import Foreign.C
import Foreign.Ptr
import Foreign.Storable

import Foreign.Marshal.Alloc
import Foreign.Marshal.Array

-- !!!!!!!!!!!!!!
-- TODO TIFFRewriteDirectory
-- !!!!!!!!!!!!!!

-- | map mxArray??
-- type Column = Int
-- type Row = Int
-- array :: (Ptr a) -> Int -> Int
array mat cols rows f = do eachColumn 0
    where
      eachColumn col | col < cols = eachRow col 0 >> eachColumn (col + 1)
                      | otherwise = return ()
      eachRow col row | row < rows = peekElemOff mat (col * rows + row) >>= \elem -> f row col elem >> eachRow col (row + 1)
                      | otherwise = return ()

array' :: (Ptr CUChar) -> Int -> Int -> (Int -> Int -> CUChar -> Int -> IO ()) -> IO ()
array' mat cols rows f = do eachColour 0
    where
      colours = 3
      eachColour colour | colour < colours = eachColumn colour 0 >> eachColour (colour + 1)
                      | otherwise = return ()
      eachColumn colour col | col < cols = eachRow colour col 0 >> eachColumn colour (col + 1)
                      | otherwise = return ()
      eachRow colour col row | row < rows = peekElemOff mat (colour * (cols * rows) + col * rows + row) 
                                          >>= 
                                            \elem -> f row col elem colour >> eachRow colour col (row + 1)
                      | otherwise = return ()

-- array2 :: (Ptr CUChar) -> Int -> Int -> (Int -> Int -> CUChar -> Int -> IO ()) -> IO ()
-- array2 mat cols rows f = array2' 0
--     where
--       s = cols * rows
--       array2' indx = let colour = indx `div` s
--                          m = indx `mod` s
--                          row = m `div` cols
--                          col = m `mod` cols
--                      in (peekElemOff mat indx) >>= \elem -> f row col elem colour

foreign export stdcall wtifc :: CString -> (Ptr CUChar) -> (Ptr CDouble) -> Int -> Int -> IO ()

type Coord = (CDouble, CDouble)


-- | TODO haddock
-- | 
wtifc :: CString -> (Ptr CUChar) -> (Ptr CDouble) -> Int -> Int -> IO ()
wtifc filename mat xy' columns' rows = do
  -- FIXME use mxGetDimensions
  let columns = columns' `div` 3
      size = (3 + 1) * rows * columns -- with alpha
  allocaBytes size $ \buf -> do
    x:y:_ <- peekArray 2 xy'
    -- TODO write mask separately
    array' mat columns rows 
               $ \i j sample colour 
                   -> let indx = (i * columns + j) 
                      in pokeElemOff buf (indx * 4 + colour) sample 
                       >> 
                         pokeElemOff buf (indx * 4 + 3) 255 -- ????
    saveTiff filename buf (x, y) columns rows

-- | Rewrite to list processing ..
printColumn :: (Ptr CUChar) -> Int -> IO ()
printColumn buf len = printColumn' buf len 0
    where 
      printColumn' buf len indx | indx < len = do { elem <- peekElemOff buf indx; printf $ {- (show indx) ++ "-nth elem: " ++ -} (show elem) ++ ", "; printColumn' buf len (indx + 1) }
                                          | otherwise = printf "\n"

  

-- | FIXME fold

fill :: (Ptr CUChar) -> Int -> IO ()
fill buf len = fill' buf len 0
               where 
                 fill' buf len indx | indx < len = do { pokeElemOff buf indx fillColor; fill' buf len (indx + 1) }
                                    | otherwise = return ()
                 fillColor = 120

saveTiff :: CString -> (Ptr CUChar) -> Coord -> Int -> Int -> IO ()
saveTiff filename buf (x, y) width length = do
  withCString "w" $ \m -> do
    image <- tIFFOpen filename m
    setFields' image
    tIFFWriteEncodedStrip image 0 (castPtr buf) (fromIntegral ((3 + 1) * width * length)) -- with alpha
    tIFFClose image
        where
          setFields' image = 
              do
                -- number of columns of pixels X number of rows of pixels in the image
                setField image tIFFTAG_IMAGEWIDTH (fromIntegral width)
                setField image tIFFTAG_IMAGELENGTH (fromIntegral length)

                setField image tIFFTAG_BITSPERSAMPLE 8 -- gray scale
                setField image tIFFTAG_SAMPLESPERPIXEL (3 + 1)

                setField image tIFFTAG_ROWSPERSTRIP (fromIntegral length)
                setField image tIFFTAG_COMPRESSION cOMPRESSION_NONE    -- Dump

                setField image tIFFTAG_PHOTOMETRIC pHOTOMETRIC_RGB

                -- setField image tIFFTAG_FILLORDER fILLORDER_MSB2LSB
                setField image tIFFTAG_PLANARCONFIG pLANARCONFIG_CONFIG
                setField image tIFFTAG_EXTRASAMPLES [eXTRASAMPLE_UNASSALPHA]

                setField image tIFFTAG_XRESOLUTION 150.0
                setField image tIFFTAG_YRESOLUTION 150.0
                setField image tIFFTAG_RESOLUTIONUNIT rESUNIT_INCH
-- ...
                setField image tIFFTAG_XPOSITION x
                printf $ "Position: " ++ (show x) ++ "x" ++ (show y) ++ "\n"
                setField image tIFFTAG_YPOSITION y


-- Local Variables: 
-- compile-command:"\"c:/Program Files/Microsoft Visual Studio 9.0/VC/BIN/nmake\""
-- End:


