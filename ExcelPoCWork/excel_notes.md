Notes on using Excel files as Data source
=========================================
#### Using “Microsoft Access Database Engine 2010 Redistributable
##### 1. ODBC driver built on OLEDB provider
- Driver name must be known exactly
- Single connection string for (*.xls, *.xlsx, *.xlsm, *.xlsb)
- Sheet names seem to need to be known precisely because they seem difficult to discover - possibly get grouped with ‘system tables’
##### 2. OLEDB provider directly
- ACE 12.0
- Different connection strings (Excel 12.0 Xml) for (*.xlsx) versus (Excel 8.0) for (*.xls) files. [Each additional extension requires a different property string]
#### Using Open XML 
- Only works with (*.xlsx) files since (*.xls) are binary not xml based
ii. *.xlsx files are zipped xml files following “Excel (.xlsx) Extensions to the Office Open XML SpreadsheetML File Format” specifications
#### Using Excel
- Build Add-in 
- Excel version dependent
- Uses Microsoft.Office.Interop.Excel
- Native excel data model (i.e. analogous to VBA calls)
