declare @count int
declare @letter_int int
declare @letter varchar(1)
declare @name varchar(12)
declare @RTS int
declare @server int
declare @top_asset int
declare @current_parent table (ID int)
declare @parent int
declare @archive int

set @RTS=1051
set @server=0
set @top_asset=1
set @letter_int=64
set @archive=0

while char(@letter_int) != 'Z'
begin
	set @letter_int = @letter_int +1
	set @count=0
	set @letter = char(@letter_int)
	set @name = 'AMB 5 ' + @letter
 
 	INSERT INTO Assets (ParentAssetID, Description, DefaultRtsID, DefaultServerID, DefaultArchiveRtsID)
 	OUTPUT inserted.AssetID INTO @current_parent
 	VALUES (@top_asset, @name, @RTS, @server, @archive)
 
 	select @parent=ID FROM @current_parent
 
 	while @count <10 
  	begin
   		set @name = 'AMB 5 ' + @letter + CAST(@count as varchar(1))
   
   		INSERT INTO Assets (ParentAssetID, Description, DefaultRtsID, DefaultServerID, DefaultArchiveRtsID)
   		VALUES (@parent, @name, @RTS, @server, @archive)
   
   		set @count = @count + 1
  	end
  
end



