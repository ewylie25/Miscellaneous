using System;
using System.Windows;

namespace WebTestUtility
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow
    {

        private readonly WebServerConfigViewModel _vm;

        public MainWindow()
        {
            Log.Information("Form - Main Window Loading");
            this.InitializeComponent();
            this._vm = new WebServerConfigViewModel();
            this.DataContext = this._vm;
        }

        private async void btnSend_Click(object sender, RoutedEventArgs e)
        {
            Log.Information("Form - Send Request Button Clicked");
            try
            {
                Log.Information("Form - Retrieving Password from form to model");
                this._vm.Password = this.pbxPassword.Password;
                Log.Information("Form - Validating password field in model");
                this._vm.Validate("Password");

                Log.Information("Form - Creating web request");
                var request = UtilityWebRequestFactory.CreateRequest(this._vm);
                Log.Information("Form - Executing Request");
                var resp = await request.ExecuteRequest();
                Log.Information("Form - Writing response to UI");
                this.WriteData(resp);
            }
            catch (Exception ex)
            {
                Log.Exception(ex);
                Log.Information("Form - Writing error to UI");
                this.WriteData(ex.ToString());
            }
        
        }


        private void WriteData(string data)
        {
            Log.Information("Form - Writing text to text block");
            this.resultsTextBlock.Text = data;
            this.groupBox1.Visibility = Visibility.Visible;
        }

        private void btnExit_Click(object sender, RoutedEventArgs e)
        {
            Log.Information("Form - Button to exit app clicked - exiting");
            this.Close();
        }

        private void BtnClearResp_OnClick(object sender, RoutedEventArgs e)
        {
            Log.Information("Form - Clearing text in text block");
            this.resultsTextBlock.Text = string.Empty;
        }
    }
}
