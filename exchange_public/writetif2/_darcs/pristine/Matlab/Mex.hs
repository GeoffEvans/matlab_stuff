{-# LANGUAGE ForeignFunctionInterface #-}
module Matlab.Mex -- (printf)
where

import Foreign.C
import Foreign.Ptr
import Foreign.Storable

type MxArray = Ptr ()

map = 123

foreign import ccall unsafe "mexPrintf" mexPrintf :: CString -> IO ()

printf :: String -> IO ()
printf xs = do withCString xs $ \s -> mexPrintf s

foreign import ccall unsafe "mxGetData" getData :: MxArray -> IO (Ptr ())

-- number of rows
foreign import ccall unsafe "mxGetM" getRows :: MxArray -> IO CUInt

-- number of columns
foreign import ccall unsafe "mxGetN" getColumns :: MxArray -> IO CUInt





