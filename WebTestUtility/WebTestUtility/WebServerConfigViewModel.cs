using System;
using System.ComponentModel;
using System.IO;

namespace WebTestUtility
{
    public enum UserOption
    {
        UseCurrentAccount = 0,
        UseSpecifiedAccount = 1
    }

    public class WebServerConfigViewModel : INotifyPropertyChanged, IDataErrorInfo
    {
        private string _requestUrl;
        private bool _useHttps;
        private string _certUrl;
        private string _certPath;
        private UserOption _userAccountOption;
        private string _userName;
        private string _password;

        public string RequestUrl
        {
            get { return this._requestUrl; }
            set
            {
                if (this._requestUrl != value)
                {
                    this._requestUrl = value;
                    this.OnPropertyChanged("RequestUrl");
                }
            }
        }

        public bool UseHttps
        {
            get { return this._useHttps; }
            set
            {
                if (this._useHttps != value)
                {
                    this._useHttps = value;
                    this.OnPropertyChanged("UseHTTPS");
                }
            }
        }

        public string CertFilePath
        {
            get { return this._certPath; }
            set
            {
                if (this._certPath != value)
                {
                    this._certPath = value;
                    this.OnPropertyChanged("CertFilePath");
                }
            }
        }

        public UserOption UserAccountOption
        {
            get { return this._userAccountOption; }
            set
            {
                if (this._userAccountOption != value)
                {
                    this._userAccountOption = value;
                    this.OnPropertyChanged("UserAccountOption");
                }
            }
        }

        public string UserName
        {
            get { return this._userName; }
            set
            {
                if (this._userName != value)
                {
                    this._userName = value;
                    this.OnPropertyChanged("UserName");
                }
            }
        }

        //Note - not binding to this property
        public string Password
        {
            get { return this._password; }
            set
            {
                if (this._password != value)
                {
                    this._password = value;
                    this.OnPropertyChanged("Password");
                }
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected void OnPropertyChanged(string propertyName)
        {
            this.PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }

        public string this[string columnName] => this.Validate(columnName);

        public string Validate(string propertyName)
        {
            // Return error message if there is error on else return empty or null string
            string validationMessage = string.Empty;
            switch (propertyName)
            {
                case "RequestUrl":
                    validationMessage = this.ValidateUrl(this.RequestUrl);
                    break;

                case "CertFilePath":
                    if (this.UseHttps)
                    {
                        validationMessage = this.ValidateFilePath(this.CertFilePath);
                    }
                    break;

                case "UserName":
                    if (this.UserAccountOption == UserOption.UseSpecifiedAccount)
                    {
                        if (string.IsNullOrWhiteSpace(this.UserName))
                        {
                            return "User Name must be specified";
                        }
                    }
                    break;
                //note - only called manually
                case "Password":
                    if (this.UserAccountOption == UserOption.UseSpecifiedAccount)
                    {
                        if (this.Password.Length == 0)
                        {
                            return "Password must be specified";
                        }
                    }
                    break;
            }
            if (!string.IsNullOrEmpty(validationMessage))
            {
                Log.Error($"Validation error: {validationMessage}");
            }
            return validationMessage;
        }

        private string ValidateFilePath(string certFilePath)
        {
            try
            {
                Uri tempUri;

                if (!Uri.TryCreate(certFilePath, UriKind.RelativeOrAbsolute, out tempUri))
                {
                    return "failed creating Uri";
                }
                if (tempUri.Scheme != Uri.UriSchemeFile)
                {
                    return "Must be File path";
                }
            }
            catch (Exception ex)
            {
                Log.Exception(ex);
                return $"Uri validation failed for unknown reasons: {ex}";
            }

            try
            {
                FileAttributes attr = File.GetAttributes(certFilePath);
                if (attr.HasFlag(FileAttributes.Directory))
                    return "Path must be a file, not directory";
            }
            catch (Exception ex)
            {
                Log.Exception(ex);
                return $"File manipulation failed: {ex}";
            }

            try
            {
                File.Open(certFilePath, FileMode.Open, FileAccess.Read).Dispose();
                return null;
            }
            catch (IOException ex)
            {
                Log.Exception(ex);
                return $"File Open failed: {ex}";
            }
            catch (Exception ex)
            {
                Log.Exception(ex);
                return $"File Open failed for unknown reasons: {ex}";
            }
        }

        private string ValidateUrl(string url)
        {
            Uri tempUri;

            if (!Uri.TryCreate(url, UriKind.Absolute, out tempUri))
            {
                return "failed creating Uri";
            }
            if (!tempUri.IsWellFormedOriginalString())
            {
                return "Malformed Uri string";
            }
            if (this.UseHttps & tempUri.Scheme != Uri.UriSchemeHttps)
            {
                return "Must be HTTPS URL";
            }
            if (!this.UseHttps & tempUri.Scheme != Uri.UriSchemeHttp)
            {
                return "Must be HTTP URL";
            }

            return null;
        }

        public string Error => "...this seems useless...";
    }
}